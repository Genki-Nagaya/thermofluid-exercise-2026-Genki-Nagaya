using Test

if !isdefined(Main, :F04NumericalDifferentiation)
    include(joinpath(@__DIR__, "run.jl"))
end
using .F04NumericalDifferentiation

@testset "F03-F04 必須テスト（配布済み）" begin
    quadratic(x) = x^2
    @test forward_difference(quadratic, 2.0, 0.5) ≈ 4.5
    @test backward_difference(quadratic, 2.0, 0.5) ≈ 3.5
    @test centered_difference(quadratic, 2.0, 0.5) ≈ 4.0

    cubic(x) = x^3
    cubic_derivative(x) = 3x^2
    study = convergence_study(cubic, cubic_derivative, 1.0, [0.2, 0.1, 0.05])
    @test all(ratio -> 1.8 <= ratio <= 2.2, study.forward_ratios)
    @test all(ratio -> 1.8 <= ratio <= 2.2, study.backward_ratios)
    @test all(ratio -> 3.9 <= ratio <= 4.1, study.centered_ratios)
end

@testset "F03-F04 自作テスト" begin
    # TODO(自作): 別の関数、評価点、または入力条件を選び、どの実装ミスを検出するか説明できるテストを一つ書く。F03単独の別テストは作らない。
    @test false
end
