# N01のテスト枠。必須テストと自作テストの両方を受講生が実装します。
# 各TODOで入力と独立した期待値を選び、配列全体・境界・安定性を確認します。

using Test

if !isdefined(Main, :N01LinearAdvection)
    include(joinpath(@__DIR__, "run.jl"))
end

@testset "N01 必須テスト" begin
    M = N01LinearAdvection

    # (1) 矩形初期条件：区間端0.5と1.0がちょうど格子点にあり，両端を含むことも確かめる
    x = [0.0, 0.5, 0.75, 1.0, 1.5]
    @test M.rectangular_initial_condition(x) == [1.0, 2.0, 2.0, 2.0, 1.0]

    # 共通の小さな入力：c=1, dt=0.5, dx=1 → C=0.5
    u_old = [1.0, 1.0, 2.0, 2.0, 1.0]
    original = copy(u_old)

    # (2) 風上：u_i - 0.5(u_i - u_{i-1})．内部 i=2,3,4 → 1.0, 1.5, 2.0．両端は旧値のまま
    up = M.upwind_step!(similar(u_old), u_old, 1.0, 0.5, 1.0)
    @test up == [1.0, 1.0, 1.5, 2.0, 1.0]
    @test u_old == original

    # (3) 中心：u_i - 0.25(u_{i+1} - u_{i-1}) → 0.75, 1.75, 2.25（風上と違い0.75と2.25が範囲外に出る）
    ce = M.centered_step!(similar(u_old), u_old, 1.0, 0.5, 1.0)
    @test ce == [1.0, 0.75, 1.75, 2.25, 1.0]
    @test u_old == original

    # (4) 境界：左端は1.0に固定，右端は右から2番目の値をコピー（勾配ゼロ）
    @test M.apply_boundary!([5.0, 2.0, 3.0, 4.0, 9.0]) == [1.0, 2.0, 3.0, 4.0, 4.0]

    # (5) 標準条件（nx=81, CFL=0.5, t=0.5）で安定・不安定を比較
    upwind = M.simulate(; scheme=:upwind)
    centered = M.simulate(; scheme=:centered)
    @test upwind.steps == centered.steps == 40
    @test upwind.dt ≈ centered.dt ≈ 0.0125
    tol = 1e-12
    @test 1.0 - tol <= upwind.minimum && upwind.maximum <= 2.0 + tol   # 初期値の範囲[1,2]に収まる
    @test centered.maximum > 2.5 && centered.minimum < 0.5             # 範囲を大きく外れて振動する
end

@testset "N01 自作テスト" begin
    # CFL=1の風上差分は u_i^{n+1} = u_{i-1}^n となり，矩形がちょうど1格子ずつずれる（数値拡散なし）．
    # t=0.5で矩形は c*t = 0.5 だけ右へ移動し，[1.0, 1.5]にあるはず
    M = N01LinearAdvection
    r = M.simulate(; scheme=:upwind, cfl=1.0)
    @test r.steps == 20
    @test r.cfl ≈ 1.0
    expected = M.rectangular_initial_condition(r.x; plateau_start=1.0, plateau_end=1.5)
    @test isapprox(r.u, expected; atol=1e-12)
end
