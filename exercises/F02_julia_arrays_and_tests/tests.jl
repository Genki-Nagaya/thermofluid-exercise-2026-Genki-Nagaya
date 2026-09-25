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
    M=F02JuliaArraysAndTests
    ints=[10, 20, 30]
    @test M.mean_temperature(ints)===20.0
    @test M.temperature_anomaly(ints)==[-10.0, 0.0, 10.0]
    @test eltype(M.temperature_anomaly(ints))==Float64

    f32=Float32[1, 2, 4]
    @test M.mean_temperature(f32) isa Float32
    @test eltype(M.temperature_anomaly(f32))==Float32

    @test_throws ArgumentError M.mean_temperature([1.0, NaN])
    @test_throws ArgumentError M.mean_temperature([1.0, Inf])
end
