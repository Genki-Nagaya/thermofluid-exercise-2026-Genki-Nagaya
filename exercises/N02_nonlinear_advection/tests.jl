# 必須テストは提供済みの期待値で確認します。末尾の自作テストだけを受講生が編集します。
# 入力の準備 → 独立した期待値との比較 → 入力保持・不正入力の確認の順に読みます。

using Test
if !isdefined(Main, :N02NonlinearAdvection)
    include(joinpath(@__DIR__, "run.jl"))
end
using .N02NonlinearAdvection
@testset "N02 必須テスト" begin
    @testset "流束と周期添字" begin
        @test burgers_flux(1.) == 0.5
        @test burgers_flux(2.) == 2.
        @test periodic_left_index(1, 4) == 4
        @test periodic_left_index(3, 4) == 2
    end
    @testset "保存形1ステップと旧配列不変" begin
        old = [1., 2., 1., 1.]
        saved = copy(old)
        new = similar(old)
        nonlinear_upwind_step!(new, old, 0.25, 1.; boundary=:periodic)
        @test new == [1., 1.625, 1.375, 1.]
        @test old == saved
    end
    @testset "定数場と境界条件" begin
        for (boundary, value) in ((:fixed, 1.), (:periodic, 0.7))
            old = fill(value, 5)
            new = similar(old)
            nonlinear_upwind_step!(new, old, 0.2, 1.; boundary)
            apply_boundary!(new; boundary)
            @test new ≈ old
        end
        @test apply_boundary!([4., 2., 3., 8.]) == [1., 2., 3., 3.]
    end
    @testset "CFL・有界性・非線形な波形変形" begin
        for (boundary, nx) in ((:fixed, 81), (:periodic, 80))
            r = simulate(; boundary, nx)
            @test r.max_cfl <= 0.5 + 1e-14
            @test 1 - 1e-13 <= minimum(r.u) <= maximum(r.u) <= 2 + 1e-13
            # c=1で同じ最終時刻まで移した矩形との離散L1差。周期では左右をつなぐ。
            shifted = [
                0.5 <=
                    (boundary == :periodic ? mod(x - r.t_final, 2.) : x - r.t_final) <=
                    1.0 ? 2. : 1. for x in r.x
            ]
            @test r.dx * sum(abs.(r.u - shifted)) > 0.2
            @test maximum(r.u[findall(>(1.5), r.x)]) > 1.2
        end
    end
end
@testset "N02 自作テスト" begin
    @testset "周期総和保存" begin
        # 非対称・非定数・非負．max|u|=2, dt/dx=0.4 → CFL=0.8 ≤ 1
        u0 = [1.0, 1.5, 2.0, 1.2, 1.0, 1.8, 1.1]
        dt, dx = 0.04, 0.1
        nsteps = 50
        old = copy(u0)
        new = similar(old)
        for _ in 1:nsteps
            nonlinear_upwind_step!(new, old, dt, dx; boundary=:periodic)
            apply_boundary!(new; boundary=:periodic)
            old, new = new, old
        end
        # 独立した期待値：初期総和（手計算で 9.6）
        @test sum(u0) ≈ 9.6
        # 丸め誤差の上限目安：1点1ステップあたり数ulp × 点数 × ステップ数
        tol = 10 * nsteps * length(u0) * eps(maximum(u0))   # ≈ 1.6e-12
        @test abs(sum(old) - sum(u0)) <= tol
        # 場が実際に変化したこと（何もしない実装で通らないように）
        @test maximum(abs.(old .- u0)) > 0.1
    end

    @testset "循環移動と更新の整合" begin
        u = [1.0, 1.8, 1.2, 2.0, 1.4]          # 非対称・非定数
        dt, dx = 0.04, 0.1                      # CFL=0.8
        S(v) = nonlinear_upwind_step!(similar(v), v, dt, dx; boundary=:periodic)
        for k in 1:(length(u)-1)              # k ≠ 0 をすべて
            @test S(circshift(u, k)) ≈ circshift(S(u), k)
        end
    end
end
