# 必須テストは、数値更新・境界条件・熱量・標準計算を順に確認します。
# 末尾の「自作テスト」では、確認したい性質と許容誤差を自分で決めます。

using Test

# 単独実行時は読み込み、Pkg.test などで読込済みなら既存のモジュールを使う。
if !isdefined(Main, :N03Diffusion)
    include(joinpath(@__DIR__, "run.jl"))
end
import .N03Diffusion

@testset "N03 必須テスト" begin
    @testset "小配列1ステップと旧配列不変" begin
        # 非対称な配列を用意し、更新前の値を別の配列に保存する。
        old = [1., 2., 4., 3.]
        saved = copy(old)
        new = similar(old)

        # 新配列の値と、旧配列が変更されていないことを確認する。
        N03Diffusion.diffusion_step!(new, old, 0.25, 1., 1.; boundary=:insulated)
        @test new == [1.5, 2.25, 3.25, 3.5]
        @test old == saved
    end

    @testset "定数場と固定端" begin
        # 各境界条件と両立する定数場は、1ステップ後も変化しない。
        for (boundary, value) in ((:insulated, 0.7), (:fixed, 0.))
            old = fill(value, 5)
            new = similar(old)
            N03Diffusion.diffusion_step!(new, old, 0.2, 1., 1.; boundary)
            @test new ≈ old
        end

        # 既定の固定温度条件では、更新後の両端が0になる。
        old = [1., 2., 4., 3.]
        new = similar(old)
        N03Diffusion.diffusion_step!(new, old, 0.25, 1., 1.)
        @test new == [0., 2.25, 3.25, 0.]
    end

    @testset "台形則の熱量" begin
        # 断熱条件の1ステップ更新前後で、端点の重みを含む熱量を比較する。
        @test N03Diffusion.thermal_content([1., 2., 4., 3.], 1.) == 8.
        @test N03Diffusion.thermal_content([1.5, 2.25, 3.25, 3.5], 1.) == 8.
    end

    @testset "標準計算の安定性" begin
        for boundary in (:fixed, :insulated)
            # 各境界条件で、既定の終了時刻まで計算する。
            r = N03Diffusion.simulate(; boundary)

            # 終了時刻・安定条件・有限性・温度範囲を確認する。
            @test r.steps * r.dt ≈ r.t_final == 1.
            @test 0 < r.fo <= 0.5
            @test all(isfinite, r.u)
            @test -1e-13 <= r.minimum <= r.maximum <= 1 + 1e-13
        end
    end
end

@testset "N03 自作テスト" begin
    @testset "断熱の台形則熱量保存（複数ステップ）" begin
        # 端点同士・端点と内部で値が異なる非対称配列（左右反転で一致しない）
        u0 = [0.3, 1.7, 0.2, 2.9, 1.1, 0.6]
        dx, D = 0.1, 1.0
        dt = 0.4 * dx^2 / D          # r = 0.4（安定範囲内）
        nsteps = 50

        # 初期熱量を run.jl とは別の書き方（区間ごとの台形の和）で計算する
        H0_indep = dx * sum((u0[i] + u0[i+1]) / 2 for i in 1:(length(u0)-1))
        @test H0_indep ≈ 0.635                       # 手計算: 0.1*(0.45+5.9)
        @test N03Diffusion.thermal_content(u0, dx) ≈ H0_indep

        # 新旧バッファを交換しながら複数ステップ進める
        u_old, u_new = copy(u0), similar(u0)
        for _ in 1:nsteps
            N03Diffusion.diffusion_step!(u_new, u_old, dt, dx, D; boundary=:insulated)
            u_old, u_new = u_new, u_old
        end

        # 許容誤差：1ステップ・1点あたり丸め eps 程度 × 点数 × ステップ数 × 値の大きさ × 領域長
        tol = nsteps * length(u0) * eps(Float64) * maximum(abs, u0) * (length(u0) - 1) * dx
        @test abs(N03Diffusion.thermal_content(u_old, dx) - H0_indep) <= tol   # tol ≈ 9.7e-14

        # 分布が実際に変化したこと（何も更新しない誤実装でも保存は成り立ってしまうため）
        @test maximum(abs, u_old .- u0) > 0.1
    end

    @testset "固定温度の正弦解析解と格子収束" begin
        D, t_final = 0.1, 1.0
        errors = Float64[]
        for nx in (41, 81, 161)
            r = N03Diffusion.simulate(; boundary=:fixed, initial=:mode, nx, fo=0.4, t_final)
            @test r.fo ≈ 0.4             # 3格子で同じ実効Fourier数（Δt ∝ Δx²）
            # 解析解を provided の関数に頼らず独立に書く
            exact = exp(-D * (pi / 2)^2 * t_final) .* sin.(pi .* r.x ./ 2)
            push!(errors, maximum(abs.(r.u .- exact)))
        end

        # 格子を細かくすると誤差が単調に減る
        @test errors[1] > errors[2] > errors[3]
        # 観測次数 p = log2(E_coarse/E_fine) が約2（実測 2.001, 2.000）
        p = [log2(errors[i] / errors[i+1]) for i in 1:2]
        @test all(1.8 .< p .< 2.2)
        # 絶対値の確認：81点で約 3.47e-5（授業ページの値）
        @test 3.0e-5 < errors[2] < 4.0e-5
    end
end
