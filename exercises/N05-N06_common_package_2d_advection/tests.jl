# N05/N06: 提供済みの必須テスト → 受講生が編集する自作テスト。必須の期待値は保持します。

using Test, ThermofluidExercise
@testset "N05 共通API" begin
    @test ThermofluidExercise.periodic_left_index(1, 4) == 4
    @test ThermofluidExercise.periodic_right_index(4, 4) == 1
    @test ThermofluidExercise.linear_flux(2., 3.) == 6.
    @test ThermofluidExercise.burgers_flux(-2.) == 2.
    @test ThermofluidExercise.fit_timestep(0.3, 1.) == (steps = 4, dt = 0.25)
end
@testset "N06 配布済み数値テスト" begin
    old = Float64[1 2 4 3; 5 7 6 8; 9 10 12 11]
    saved = copy(old)
    new = similar(old)
    N = ThermofluidExercise.N06
    @test N.advection_step!(new, old, 0.1, 0.5, 0.25; cx = 1., cy = 0.75)===new
    @test new≈[3.2 3.3 5.0 4.9; 5.1 5.4 5.9 6.4; 8.8 9.1 10.2 10.7] atol = 1e-13
    @test old == saved
    @test N.stable_timestep(1., 0.75, 0.5, 0.25)≈0.16
    @test_throws ArgumentError N.advection_step!(old, old, 0.1, 0.5, 0.25)
    @test_throws ArgumentError N.advection_step!(
        new,
        old,
        0.21,
        0.5,
        0.25;
        cx = 1.,
        cy = 0.75,
    )
    N.advection_step!(new, fill(2., 3, 4), 0.1, 0.5, 0.25)
    @test new == fill(2., 3, 4)
    N.advection_step!(new, old, 0.1, 0.5, 0.25; cx = 0., cy = 0.)
    @test new == old
    N.advection_step!(new, old, 0.1, 0.5, 0.25; cx = 0., cy = 0.75)
    @test new≈[1.6 1.7 3.4 3.3; 5.9 6.4 6.3 7.4; 9.6 9.7 11.4 11.3]
end
@testset "N06 自作: 複数ステップの保存量" begin
    # x,yで形も格子数も異なる非定数場を20ステップ進める。周期境界の保存形なので総和は不変、
    # 合成CFL<=1では更新が旧値の凸結合になるため、最小・最大は初期値の範囲を出ない。
    N = ThermofluidExercise.N06
    nx, ny, dx, dy = 7, 5, 0.3, 0.2
    u0 = [sin(0.9i) + 0.3j^2 + 0.1i * j for i in 1:nx, j in 1:ny]
    saved = copy(u0)
    dt = N.stable_timestep(1.3, 0.4, dx, dy)
    a, b = copy(u0), similar(u0)
    for _ in 1:20
        N.advection_step!(b, a, dt, dx, dy; cx = 1.3, cy = 0.4)
        a, b = b, a
    end
    @test u0 == saved
    @test a != u0
    @test sum(a)≈sum(u0) rtol = 1e-13
    @test minimum(u0) - 1e-13 <= minimum(a) && maximum(a) <= maximum(u0) + 1e-13
    # Cx=1, Cy=0 では風上更新が1格子の厳密なシフトになり、x方向の周期左隣を取り違えると一致しない。
    N.advection_step!(b, u0, dx / 1.3, dx, dy; cx = 1.3, cy = 0.)
    @test b≈circshift(u0, (1, 0)) atol = 1e-13
    N.advection_step!(b, u0, dy / 0.4, dx, dy; cx = 0., cy = 0.4)
    @test b≈circshift(u0, (0, 1)) atol = 1e-13
end
include("simulate.jl")
@testset "N06 自作: 3格子の解析解誤差と収束" begin
    # 公式3格子をt=1まで計算し、解析解との二乗平均誤差を比べる。
    # 一次風上差分なので、格子を半分にするごとに誤差はおよそ半分（観測次数≈1）になるはずである。
    S = N06Simulation
    errors = Float64[]
    for (nx, ny) in S.OFFICIAL_GRIDS
        c = S.simulate_case(nx, ny)
        @test c["time"][end] == 1.
        exact = S.exact_field(c["x"], c["y"], 1.)
        push!(errors, sqrt(sum(abs2, c["u"][:, :, end] - exact) / (nx * ny)))
    end
    orders = log2.(errors[1:(end - 1)] ./ errors[2:end])
    @test all(>(0), errors)
    @test all(diff(errors) .< 0)
    @test all(0.8 .<= orders .<= 1.2)
    println("N06 3格子: L2誤差 = ", errors, ", 観測次数 = ", orders)
end
module N0506OutputChecks
include("provided_support.jl")
include("N05.jl")
end
@testset "N05-N06 保存結果と出自" begin
    @test N0506OutputChecks.check_complete()
    baseline = joinpath(@__DIR__, "results", "N05", "baseline.toml")
    regression = joinpath(@__DIR__, "results", "N05", "regression.toml")
    @test isfile(baseline) && isfile(regression)
    report = N0506OutputChecks.TOML.parsefile(regression)
    @test report["baseline_sha256"] == N0506OutputChecks.file_sha(baseline)
    @test report["code_sha256"] == N0506OutputChecks.N05Regression.code_hashes(
        normpath(joinpath(@__DIR__, "..", "..")),
    )
end
