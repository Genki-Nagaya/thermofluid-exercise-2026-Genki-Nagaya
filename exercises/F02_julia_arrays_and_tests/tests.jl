using Test

if !isdefined(Main, :F02JuliaArraysAndTests)
    include(joinpath(@__DIR__, "run.jl"))
end

@testset "F02 必須テスト（配布済み）" begin
    values = [5.0, 7.0, 12.0]
    original = copy(values)
    anomalies = F02JuliaArraysAndTests.temperature_anomaly(values)

    @test F02JuliaArraysAndTests.mean_temperature(values) == 8.0
    @test anomalies == [-3.0, -1.0, 4.0]
    @test isapprox(sum(anomalies), 0.0; atol=100eps())
    @test values == original
    @test_throws ArgumentError F02JuliaArraysAndTests.mean_temperature(Float64[])
end

@testset "F02 自作テスト" begin
    # TODO(自作): 戻り値の型、別の数学的性質、または必須とは異なる不正入力から一つ選び、入力と期待値を自分で書く。
    @test false
end
