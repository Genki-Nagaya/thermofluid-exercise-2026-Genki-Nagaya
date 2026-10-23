# 必須テストは提供済みの期待値で確認します。末尾の自作テストだけを受講生が編集します。
# 入力の準備 → 独立した期待値との比較 → 入力保持・不正入力の確認の順に読みます。

using Test
if !isdefined(Main, :N04AdvectionDiffusion)
    include(joinpath(@__DIR__, "run.jl"))
end
import .N04AdvectionDiffusion
const N04 = N04AdvectionDiffusion
N04.validate_model(N04.SELECTED_MODEL)
@testset "N04 必須テスト" begin
    model = N04.SELECTED_MODEL
    @testset "非対称1ステップ・周期端・旧配列不変" begin
        old = [1., 2., 4., 3.]
        saved = copy(old)
        new = similar(old)
        expected = model == :linear ? [1.26, 1.92, 3.74, 3.08] : [1.46, 1.87, 3.34, 3.33]
        @test N04.advection_diffusion_step!(new, old, 0.1, 1., 0.2) === new
        @test new ≈ expected atol = 1e-14
        @test old == saved
    end
    @testset "定数場と周期積分" begin
        old = fill(1.7, 5)
        new = similar(old)
        N04.advection_diffusion_step!(new, old, 0.1, 1., 0.2)
        @test new ≈ old
        @test N04.conserved_integral([1., 2., 4., 3.], 0.5) == 5.
    end
    @testset "移流・拡散をゼロにした極限" begin
        old = [1., 2., 4., 3.]
        new = similar(old)
        N04.advection_diffusion_step!(new, old, 0.1, 1., 0.2; advection=false)
        @test new ≈ [1.06, 2.02, 3.94, 2.98]
        N04.advection_diffusion_step!(new, old, 0.1, 1., 0.)
        @test new ≈ (model == :linear ? [1.2, 1.9, 3.8, 3.1] : [1.4, 1.85, 3.4, 3.35])
    end
    @testset "合成刻み・最終時刻・標準安定性" begin
        @test N04.stable_timestep(2., 0.5, 0.1) ≈ 1 / 6
        r = N04.simulate()
        @test r.steps * r.dt ≈ r.t_final == 1.
        @test r.max_stability_number <= 0.8 + 1e-12
        @test all(isfinite, r.u) && 1 - 1e-12 <= r.minimum <= r.maximum <= 2 + 1e-12
    end
end
@testset "N04 自作テスト" begin
    model = N04.SELECTED_MODEL

    @testset "複数ステップの周期積分保存" begin
        # 非負・非定数・非対称（周期シフトでも左右反転でも一致しない）
        u0 = [1.2, 1.9, 1.05, 1.7, 1.4, 1.1]
        dx, D, dt, nsteps = 0.25, 0.1, 0.05, 20
        # 合成安定数：線形 C+2r=0.2+0.16=0.36，非線形 0.38+0.16=0.54（どちらも≤1）

        # 初期積分を手計算で独立に与える：0.25 × 8.35
        I0 = 2.0875
        @test N04.conserved_integral(u0, dx) ≈ I0

        u_old, u_new = copy(u0), similar(u0)
        for _ in 1:nsteps
            N04.advection_diffusion_step!(u_new, u_old, dt, dx, D; model)
            u_old, u_new = u_new, u_old
        end

        # 許容誤差：ステップ数 × 点数 × eps × 値の大きさ × 領域長 ≈ 7.6e-14
        tol = nsteps * length(u0) * eps(Float64) * maximum(abs, u0) * length(u0) * dx
        @test abs(dx * sum(u_old) - I0) <= tol

        # 何もしない実装でも保存は成り立つので，分布が実際に変化したことも確認する
        @test maximum(abs, u_old .- u0) > 0.1
        # 最大値原理：初期の範囲 [1.05, 1.9] に収まる
        @test 1.05 - 1e-12 <= minimum(u_old) && maximum(u_old) <= 1.9 + 1e-12
    end

    @testset "解析解と格子収束（40・80・160点）" begin
        D, t_final, k = 0.1, 1.0, 2pi / 2
        # 解析解を provided の関数に頼らず独立に書く
        exact(x, t) = model == :linear ?
                      1 .+ 0.5exp(-D * k^2 * t) .* sin.(k .* (x .- t)) :
                      (A=0.5exp(-D * k^2 * t); θ=k .* (x .- t);
            1 .+ 2D * k * A .* sin.(θ) ./ (1 .+ A .* cos.(θ)))

        errors = Float64[]
        for nx in (40, 80, 160)
            r = N04.simulate(; model, nx, initial=:mode, t_final)
            @test r.steps * r.dt ≈ t_final          # 同じ最終時刻で比較
            push!(errors, maximum(abs.(r.u .- exact(r.x, t_final))))
        end

        @test errors[1] > errors[2] > errors[3]
        p = [log2(errors[i] / errors[i+1]) for i in 1:2]
        @test all(0.8 .< p .< 1.2)                  # 課題の受入範囲
        @test p[2] > p[1]                           # 細かくするほど1次へ近づく
        # 絶対値：80点の誤差（線形 1.994e-2，非線形 1.340e-2）の±10%以内
        ref80 = model == :linear ? 1.994e-2 : 1.340e-2
        @test 0.9 * ref80 < errors[2] < 1.1 * ref80
    end
end
