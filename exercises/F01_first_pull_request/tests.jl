using Test

if !isdefined(Main, :F01FirstPullRequest)
    include(joinpath(@__DIR__, "run.jl"))
end

@testset "F01 必須テスト（配布済み）" begin
    @test F01FirstPullRequest.student_greeting("  Thermofluid  ") ==
        "Hello, Thermofluid!"
    @test_throws ArgumentError F01FirstPullRequest.student_greeting("   ")
end

@testset "F01 自作テスト" begin
    "Hello, Euler!"
end
