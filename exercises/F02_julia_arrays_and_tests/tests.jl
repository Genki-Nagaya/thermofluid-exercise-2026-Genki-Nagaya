# F02の課題テスト。提供済みの必須確認 → 受講生が追加する自作テストの順に読みます。
# 必須テストの期待値は保持し、自作テストのTODOを編集します。

using Test

if !isdefined(Main, :F02JuliaArraysAndTests)
    include(joinpath(@__DIR__, "run.jl"))
end

# 提供済みの期待値を使い、実装と独立した基準で確認する。
@testset "F02 必須テスト（配布済み）" begin
    values = [5.0, 7.0, 12.0]
    original = copy(values)
    anomalies = F02JuliaArraysAndTests.temperature_anomaly(values)

    # この具体例の平均と偏差は厳密に表せる値なので、==で比較する。
    @test F02JuliaArraysAndTests.mean_temperature(values) == 8.0
    @test anomalies == [-3.0, -1.0, 4.0]

    # originalは呼び出し前にcopyした値。単なる代入では変更を見逃す。
    @test values == original
    @test_throws ArgumentError F02JuliaArraysAndTests.mean_temperature(Float64[])
end

# 必須テストが扱わない条件を選び、期待値を自分で決める。
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
