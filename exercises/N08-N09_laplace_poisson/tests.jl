# N08/N09: 提供済みの必須確認 → 受講生が編集する自作テスト。

using Test
# 実行対象はN08-checkまたはall。
const N08N09_TEST_MODE =
    abspath(PROGRAM_FILE) == (@__FILE__) ?
    (isempty(ARGS) ? "all" : length(ARGS) == 1 ? only(ARGS) : "invalid") : "all"
N08N09_TEST_MODE in ("N08-check", "all") || error("引数は[N08-check|all]です")
using ThermofluidExercise
const E = ThermofluidExercise.Elliptic
include("provided_support.jl")
@testset "N08: 小配列・境界・独立残差・残差停止" begin
    old = [Float64(i ^ 2 + 3j ^ 2 + i * j) for i in 1:5, j in 1:4]
    saved = copy(old)
    new = fill(-7., 5, 4)
    g = [Float64(10i + j) for i in 1:5, j in 1:4]
    inside = copy(new[2:4, 2:3])
    @test E.apply_dirichlet!(new, g)===new
    @test new[1, :] == g[1, :] &&
          new[end, :] == g[end, :] &&
          new[:, 1] == g[:, 1] &&
          new[:, end] == g[:, end]
    @test new[2:4, 2:3] == inside
    @test E.laplace_jacobi_step!(new, old, 0.2, 0.3)===new
    @test new[2:4, 2:3]≈[281 / 13 502 / 13; 372 / 13 606 / 13; 489 / 13 736 / 13]
    @test old == saved
    @test new[1, :] == old[1, :] &&
          new[end, :] == old[end, :] &&
          new[:, 1] == old[:, 1] &&
          new[:, end] == old[:, end]
    r = similar(old)
    @test E.laplace_residual!(r, old, 0.2, 0.3)===r
    @test r[2:4, 2:3]≈fill(-350 / 3, 3, 2)
    @test all(iszero, r[1, :]) &&
          all(iszero, r[end, :]) &&
          all(iszero, r[:, 1]) &&
          all(iszero, r[:, end])
    @test E.residual_converged(1., 10.; atol = 1., rtol = 0.01)
    @test !E.residual_converged(1.01, 10.; atol = 1., rtol = 0.01)
    @test E.residual_converged(0., 0.; atol = 1e-10, rtol = 0.)
    zero = E.solve_laplace(zeros(5, 4), zeros(5, 4), 0.2, 0.3)
    @test zero.converged && zero.iterations == 0 && zero.update_history == [0.]
    boundary = zeros(5, 4)
    boundary[end, :].=1e-12
    # 更新量と方程式残差の違いを確認する。
    one = E.solve_laplace(
        zeros(5, 4),
        boundary,
        1e-7,
        2e-7;
        maxiter = 1,
        atol = 1e-10,
        rtol = 0.,
    )
    @test one.reason == :maxiter && !one.converged && one.iterations == 1
    @test one.update_history[end] < 1e-10 && one.residual_history[end] > one.threshold
    independent = independent_residual(one.u, zeros(5, 4), 1e-7, 2e-7)
    @test independent≈one.residual_history[end]
    @test length(one.residual_history) == length(one.update_history) == 2
    @test_throws ArgumentError E.laplace_jacobi_step!(old, old, 0.2, 0.3)
end
if N08N09_TEST_MODE == "N08-check"
    println("N08授業内検証成功")
else
    @testset "N09: 生成項の符号・f=0一致・残差" begin
        old = [Float64(i ^ 2 + 3j ^ 2 + i * j) for i in 1:5, j in 1:4]
        saved = copy(old)
        f = fill(2., 5, 4)
        new = similar(old)
        r = similar(old)
        lap = similar(old)
        @test E.poisson_jacobi_step!(new, old, f, 0.2, 0.3)===new
        @test new[2:4, 2:3]≈[
            7016 / 325 12541 / 325;
            9291 / 325 15141 / 325;
            12216 / 325 18391 / 325
        ]
        @test old == saved && new[1, :] == old[1, :] && new[:, end] == old[:, end]
        @test E.poisson_residual!(r, old, f, 0.2, 0.3)===r
        @test r[2:4, 2:3]≈fill(-344 / 3, 3, 2)
        @test all(iszero, r[1, :]) && all(iszero, r[:, end])
        E.laplace_jacobi_step!(lap, old, 0.2, 0.3)
        E.poisson_jacobi_step!(new, old, zeros(5, 4), 0.2, 0.3)
        @test lap == new
        E.laplace_residual!(lap, old, 0.2, 0.3)
        E.poisson_residual!(r, old, zeros(5, 4), 0.2, 0.3)
        @test lap == r
    end
    @testset "自作1: 非対称小配列・両問題1 sweep・旧場・境界・f=0一致" begin
        # TODO: 手計算した入力と期待値を使い、配布済みassertionを保持して独自に検証する。
        # 入力: 4x3の非対称配列、dx=0.5, dy=1 → ax=4, ay=1, d=10（ax≠ayなので軸の取り違えを検出できる）。
        # 内部点は(2,2)と(3,2)の2点。(3,2)の西隣(2,2)は旧値5を使うので、Gauss-Seidel的な混入も検出できる。
        old = [
            0. 1. 0.
            2. 5. 6.
            1. 4. 2.
            0. 0. 1.
        ]
        saved = copy(old)
        # 境界のfは使われないはずなので、大きな値を置いて誤用を検出する。
        f = [
            100. 100. 100.
            100. 3. 100.
            100. -2. 100.
            100. 100. 100.
        ]
        fsaved = copy(f)
        boundary_equal(a, b) =
            a[1, :] == b[1, :] && a[end, :] == b[end, :] && a[:, 1] == b[:, 1] && a[:, end] == b[:, end]
        # Laplace 1 sweep（手計算）:
        #   (2,2): (4*(4+1) + 1*(6+2))/10 = 2.8   （軸を取り違えると3.7）
        #   (3,2): (4*(0+5) + 1*(2+1))/10 = 2.3   （更新済み2.8を使うと(4*2.8+3)/10=1.42）
        lap_new = fill(NaN, 4, 3)
        @test E.laplace_jacobi_step!(lap_new, old, 0.5, 1.0)===lap_new
        @test lap_new[2:3, 2]≈[2.8, 2.3] rtol = 1e-14
        @test boundary_equal(lap_new, old) && old == saved
        # Poisson 1 sweep: u = (… - f)/d
        #   (2,2): (28 - 3)/10 = 2.5、(3,2): (23 + 2)/10 = 2.5   （fの符号を逆にすると3.1と2.1）
        poi_new = fill(NaN, 4, 3)
        @test E.poisson_jacobi_step!(poi_new, old, f, 0.5, 1.0)===poi_new
        @test poi_new[2:3, 2]≈[2.5, 2.5] rtol = 1e-14
        @test boundary_equal(poi_new, old) && old == saved && f == fsaved
        # 残差 r = f - Δ_h u（手計算）:
        #   (2,2): Δ_h u = 4*(4-10+1) + 1*(6-10+2) = -22 → Laplace 22、Poisson 3+22 = 25
        #   (3,2): Δ_h u = 4*(0-8+5) + 1*(2-8+1) = -17 → Laplace 17、Poisson -2+17 = 15
        r_lap = fill(NaN, 4, 3)
        r_poi = fill(NaN, 4, 3)
        @test E.laplace_residual!(r_lap, old, 0.5, 1.0)===r_lap
        @test E.poisson_residual!(r_poi, old, f, 0.5, 1.0)===r_poi
        @test r_lap[2:3, 2]≈[22., 17.] rtol = 1e-14
        @test r_poi[2:3, 2]≈[25., 15.] rtol = 1e-14
        # 境界残差は0（NaNで初期化したので、書き忘れも検出する）。
        @test boundary_equal(r_lap, zeros(4, 3)) && boundary_equal(r_poi, zeros(4, 3))
        @test old == saved && f == fsaved
        # f=0ならPoissonはLaplaceに一致する。
        z = zeros(4, 3)
        poi0 = similar(old)
        r0 = similar(old)
        E.poisson_jacobi_step!(poi0, old, z, 0.5, 1.0)
        E.poisson_residual!(r0, old, z, 0.5, 1.0)
        @test poi0 == lap_new && r0 == r_lap
        # apply_dirichlet!は指定値gを四辺へ書き、内部点は保持する。
        u = copy(old)
        g = [Float64(10i + j) for i in 1:4, j in 1:3]
        E.apply_dirichlet!(u, g)
        @test boundary_equal(u, g) && u[2:3, 2] == old[2:3, 2]
    end
    @testset "自作2: 両問題の残差停止・解析解誤差・3格子の両区間収束" begin
        # TODO: LaplaceとPoissonの両方、p=log2(E_coarse/E_fine)を検証する。
        # 公式出力(17x13〜65x49)とは別の粗い3格子 9x5, 17x9, 33x17 を使う（dx=dy=0.25→0.125→0.0625）。
        # 期待値: 5点差分は二次精度なので、格子幅を半分にするとLinf誤差は約1/4、p≈2。
        # 許容: 1.8 <= p <= 2.2（粗格子の高次項と、反復停止による残りの誤差の分だけ2からずれる）。
        grids = ((9, 5), (17, 9), (33, 17))
        for id in ("N08", "N09")
            errs = Float64[]
            for (nx, ny) in grids
                p = problem(id, nx, ny)
                res =
                    id == "N08" ? E.solve_laplace(p.u0, p.g, p.dx, p.dy) :
                    E.solve_poisson(p.u0, p.g, p.f, p.dx, p.dy)
                # 残差停止: 収束で止まり、最後の残差は閾値以下、その直前は閾値より大きい（最初に下回った回で止まる）。
                @test res.converged && res.reason == :converged
                @test res.residual_history[end] <= res.threshold < res.residual_history[end - 1]
                @test res.threshold == max(1e-10, 1e-12 * res.residual_history[1])
                # 提供の独立残差で、保存される場が閾値を満たすことを確かめる。
                f = id == "N08" ? zeros(nx, ny) : p.f
                @test independent_residual(res.u, f, p.dx, p.dy) <= res.threshold
                # 境界は指定値のまま。
                @test res.u[1, :] == p.g[1, :] && res.u[:, end] == p.g[:, end]
                push!(errs, maximum(abs, res.u .- p.exact))
            end
            # 両区間（粗→中、中→細）で誤差が減り、観測次数が2に近い。
            orders = [log2(errs[k] / errs[k + 1]) for k in 1:2]
            @test errs[1] > errs[2] > errs[3]
            @test all(p->1.8 <= p <= 2.2, orders)
        end
    end
    @testset "公式12出力・来歴・容量" begin
        @test check_complete()
    end
end
