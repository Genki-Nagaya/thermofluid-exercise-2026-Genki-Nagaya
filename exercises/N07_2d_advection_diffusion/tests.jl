# N07: 提供済みの独立した期待値による必須確認 → 受講生が編集する自作テスト。

module N07Tests
using Test, ThermofluidExercise
const N = ThermofluidExercise.N07
include("provided_support.jl")
include("analyze.jl")
# 配布済みの必須テスト。自作欄と区別し、入力・期待値を保持する。
@testset "N07 independent faces, half-cell walls and heat" begin
    a = Float64[1 2 4 3; 5 7 6 8; 9 10 12 11]
    saved = copy(a)
    b = similar(a)
    bc = N.channel_boundaries(:dirichlet)
    f = N.flux_buffers(a)
    @test N.thermal_stable_timestep(1., 0., 0.05, 0.5, 0.25, bc)≈0.8 / 5.0
    @test N.thermal_fluxes!(f, a, 0.5, 0.25; cx = 1., cy = 0., kappa = 0.05, bc)===f
    @test f.adv_x == [1 1 1 1; 1 2 4 3; 5 7 6 8; 9 10 12 11]
    @test f.diff_x≈[0 -0.2 -0.6 -0.4; -0.4 -0.5 -0.2 -0.5; -0.4 -0.3 -0.6 -0.3; 0 0 0 0]
    @test f.adv_y == zeros(3, 5)
    @test f.diff_y≈[-0.4 -0.2 -0.4 0.2 1.2; -2 -0.4 0.2 -0.4 3.2; -3.6 -0.2 -0.4 0.2 4.4]
    r = N.thermal_step!(b, a, 0.01, 0.5, 0.25; cx = 1., cy = 0., kappa = 0.05, bc)
    @test r.temperature===b
    @test collect(r.boundary_rates.advective) == [1., -10.5, 0., 0.]
    @test collect(r.boundary_rates.diffusive)≈[-0.3, 0., -3., -4.4]
    @test b≈[1.0 1.994 3.908 2.922; 4.856 6.872 5.992 7.752; 8.776 9.942 11.844 10.766] atol =
        1e-14
    @test 0.5 * 0.25 * sum(b - a)≈0.01 * (-17.2) atol = 1e-14
    @test a == saved
    for wall in (:dirichlet, :insulated)
        bc = N.channel_boundaries(wall)
        q0 = 0.125sum(a)
        total = 0.
        u = copy(a)
        v = similar(u)
        for _ in 1:20
            r = N.thermal_step!(v, u, 0.01, 0.5, 0.25; cx = 1., cy = 0., kappa = 0.05, bc)
            rates = r.boundary_rates
            total+=0.01 * (sum(rates.advective) + sum(rates.diffusive))
            u, v = v, u
        end
        @test 0.125sum(u) - q0≈total atol = 1e-12
        wall == :insulated &&
            @test r.boundary_rates.diffusive.south == r.boundary_rates.diffusive.north == 0.
    end
    for kind in (:periodic, :insulated, :dirichlet)
        bc = N.boundaries(kind; value = kind == :dirichlet ? 2. : 0.)
        N.thermal_step!(
            b,
            fill(2., 3, 4),
            0.01,
            0.5,
            0.25;
            cx = 0.,
            cy = 0.,
            kappa = 0.05,
            bc,
        )
        @test b≈fill(2., 3, 4)
    end
    bc = N.boundaries(:dirichlet)
    bound = 1 / (0.05 * (3 / 0.5 ^ 2 + 3 / 0.25 ^ 2))
    impulse = zeros(3, 4)
    impulse[1, 1] = 1.
    N.thermal_step!(b, impulse, bound, 0.5, 0.25; cx = 0., cy = 0., kappa = 0.05, bc)
    @test minimum(b) >= -1e-15
    fill!(b, -99.)
    @test_throws ArgumentError N.thermal_step!(
        b,
        a,
        bound * (1 + 1e-8),
        0.5,
        0.25;
        cx = 0.,
        cy = 0.,
        kappa = 0.05,
        bc,
    )
    @test all(==(-99.), b)
    @test_throws ArgumentError N.thermal_step!(
        b,
        a,
        1 / (0.1 * (1 / 0.5 ^ 2 + 1 / 0.25 ^ 2)),
        0.5,
        0.25;
        cx = 0.,
        cy = 0.,
        kappa = 0.05,
        bc,
    )
    # Signed periodic flow at both ends, old-state use and heat cancellation.
    bc = N.boundaries(:periodic)
    for (cx, cy, kappa) in ((-1., 0.5, 0.), (0., 0., 0.05), (0., 0., 0.))
        r = N.thermal_step!(b, a, 0.01, 0.5, 0.25; cx, cy, kappa, bc)
        @test sum(b)≈sum(a) atol = 1e-12
        @test sum(r.boundary_rates.advective) + sum(r.boundary_rates.diffusive)≈0. atol =
            1e-13
        if kappa == 0 && cx != 0
            @test b[3, 1]≈8.88
        elseif cx == cy == kappa == 0
            @test b == a
        end
    end
end
@testset "N07 two-component signed advective Burgers" begin
    u = [1. -2. 0. 0.5; -0.5 2. -1. 0.; 3. -1. 0.2 -2.]
    v = [-0.5 1. 2. -1.; 2. -1. 0. 0.5; 0. 0.5 -2. 1.]
    originals = (copy(u), copy(v))
    un = similar(u)
    vn = similar(v)
    dt = 0.001
    dx = 0.5
    dy = 0.25
    nu = 0.05
    r = N.burgers_step!(un, vn, u, v, dt, dx, dy; nu)
    @test r.u===un && r.v===vn
    @test un≈[
        0.9953 -1.967 -0.01736 0.4989;
        -0.4891 1.9662 -0.99396 -0.0035;
        2.9707 -0.98944 0.17892 -1.97454
    ] atol = 1e-14
    @test vn≈[
        -0.4946 0.9851 1.9876 -0.9925;
        1.9815 -0.9841 -0.0044 0.4996;
        0.0135 0.4974 -1.9696 0.9763
    ] atol = 1e-14
    @test u == originals[1] && v == originals[2]
    @test N.burgers_stable_timestep(u, v, nu, dx, dy)≈0.8 / (9.0 + 2.0)
    @test N.burgers_stable_timestep(2u, 2v, nu, dx, dy) <
          N.burgers_stable_timestep(u, v, nu, dx, dy)
    a = fill(2., 3, 4)
    b = fill(-1., 3, 4)
    N.burgers_step!(un, vn, a, b, 0.001, dx, dy; nu = 0.)
    @test un == a && vn == b
    N.burgers_step!(un, vn, u, zeros(3, 4), dt, dx, dy; nu = 0.)
    @test all(iszero, vn)
    arrays = [copy(u) for _ in 1:4]
    for (i, j) in ((1, 3), (1, 2))
        args = copy(arrays)
        args[j] = args[i]
        @test_throws ArgumentError N.burgers_step!(args..., dt, dx, dy; nu)
    end
    fill!(un, -99.)
    fill!(vn, -88.)
    @test_throws ArgumentError N.burgers_step!(un, vn, u, v, 1., dx, dy; nu)
    @test all(==(-99.), un) && all(==(-88.), vn)
    @test_throws ArgumentError N.burgers_stable_timestep(u, v, -1., dx, dy)
end

@testset "N07 配布済み熱収支解析" begin
    b = N07Analysis.heat_budget(
        [0., 0.25, 1.],
        [2., 4., 5.],
        [3. 4.; -1. -2.; 0. 0.; 0. 0.],
        [0.5 0.2; 0. 0.; -0.2 -0.4; -0.3 -0.8],
    )
    @test b.net_input≈[2., 1.]
    @test b.cumulative_input≈[0., 2., 3.]
    @test b.residual≈zeros(3) atol = 1e-14
    @test_throws ArgumentError N07Analysis.heat_budget(
        [0., 0.],
        [1., 1.],
        zeros(4, 1),
        zeros(4, 1),
    )
end
@testset "N07 自作: 開境界2条件の辺別熱収支" begin
    # TODO(自作): 複数ステップで辺別・移流/拡散別の積算、熱量変化、断熱壁を検証する。
    # 入力と独立期待値、許容値の根拠を記す。公式出力は上書きしない。
    # 入力: 非対称な5×4場、dx≠dy、左流入T_in=1・右流出・上下壁(固定0/断熱)。15ステップ。
    # 独立期待値: 各stepの旧温度から、課題表の境界流束を流入正・辺長つきでテスト内に直接書いた式。
    # 許容値: 値はO(1)、20セル×15stepの加算の丸めは1e-14程度なので、1e-12で十分に余裕がある。
    T0 = [0.2 0.5 0.9 0.3; 0.7 0.1 0.4 0.8; 0.6 0.95 0.25 0.45; 0.35 0.75 0.55 0.15; 0.85 0.05 0.65 0.5]
    dx, dy, cx, kappa, Tin, nsteps = 0.3, 0.2, 1.2, 0.04, 1.0, 15
    nx, ny = size(T0)
    for wall in (:dirichlet, :insulated)
        bc = N.channel_boundaries(wall)
        dt = N.thermal_stable_timestep(cx, 0., kappa, dx, dy, bc)
        # 半セルの係数3が両軸に入るかどうかで刻みが変わることも確かめる。
        ay = wall == :dirichlet ? 3 : 2
        @test dt ≈ 0.8 / (cx / dx + kappa * (3 / dx^2 + ay / dy^2))
        u = copy(T0)
        v = similar(u)
        got = Dict((m, s) => 0.0 for m in (:advective, :diffusive), s in N.SIDES)
        want = copy(got)
        heat = [dx * dy * sum(u)]
        for _ in 1:nsteps
            # 表の外向き流束を内向き正へ反転し、面の長さ(西東dy・南北dx)を掛けた独立の期待値。
            want[(:advective, :west)] += dt * dy * ny * cx * Tin
            want[(:advective, :east)] += dt * (-dy * cx * sum(u[nx, :]))
            want[(:diffusive, :west)] += dt * dy * sum(2kappa * (Tin .- u[1, :]) / dx)
            if wall == :dirichlet
                want[(:diffusive, :south)] += dt * (-dx * sum(2kappa * u[:, 1] / dy))
                want[(:diffusive, :north)] += dt * (-dx * sum(2kappa * u[:, ny] / dy))
            end
            r = N.thermal_step!(v, u, dt, dx, dy; cx, cy = 0., kappa, bc)
            for m in (:advective, :diffusive), s in N.SIDES
                got[(m, s)] += dt * getproperty(getproperty(r.boundary_rates, m), s)
            end
            # 断熱壁・流出端の拡散と、上下壁の移流は毎stepちょうど0。
            @test r.boundary_rates.diffusive.east == 0.0
            @test r.boundary_rates.advective.south == r.boundary_rates.advective.north == 0.0
            wall == :insulated &&
                @test r.boundary_rates.diffusive.south == r.boundary_rates.diffusive.north == 0.0
            u, v = v, u
            push!(heat, dx * dy * sum(u))
        end
        for k in keys(want)
            @test got[k] ≈ want[k] atol = 1e-12
        end
        # 辺別・機構別の累積の総和が熱量変化と一致する（熱収支）。
        @test heat[end] - heat[1] ≈ sum(values(got)) atol = 1e-12
        # 固定温度の壁は正温度の場から熱を奪い、断熱との違いが出る。
        wall == :dirichlet && @test got[(:diffusive, :south)] < -1e-3
        # 区間積分(5step×3区間)をheat_budgetへ渡し、残差が丸め誤差内であることも確認する。
        adv = zeros(4, 3)
        dif = zeros(4, 3)
        u = copy(T0)
        for k in 1:3, _ in 1:5
            r = N.thermal_step!(v, u, dt, dx, dy; cx, cy = 0., kappa, bc)
            adv[:, k] .+= dt .* collect(r.boundary_rates.advective)
            dif[:, k] .+= dt .* collect(r.boundary_rates.diffusive)
            u, v = v, u
        end
        b = N07Analysis.heat_budget(collect(0:3) .* 5dt, heat[1:5:16], adv, dif)
        @test maximum(abs, b.residual) <= 1e-12
        @test b.cumulative_input[end] ≈ sum(values(got)) atol = 1e-12
    end
end
@testset "N07 自作: 温度とBurgersの解析解収束" begin
    # TODO(自作): 周期温度の拡散・合成と二成分Burgersを3格子で計算し、誤差と両区間の次数を検証する。
    # 公式出力を使わず、t=0.5(公式は1)までテスト内で時間ループを回す。誤差はRMS(課題のE_T,E_B)。
    # 期待値: 拡散と閉じた固定温度は2次(1.8<=p<=2.2)、移流を含む温度とBurgersは1次風上なので0.8<=p<=1.2。
    grids = ((40, 30), (80, 60), (160, 120))
    tf = 0.5
    function temperature_error(id, nx, ny)
        c = thermal_config(id)
        dx, dy = 2 / nx, 1 / ny
        x, y = cell_centers(nx, 2.), cell_centers(ny, 1.)
        a = exact_temperature(x, y, 0., id)
        b = similar(a)
        dtmax = N.thermal_stable_timestep(c.cx, c.cy, c.kappa, dx, dy, c.bc)
        t = 0.
        while t < tf
            dt = min(dtmax, tf - t)
            N.thermal_step!(b, a, dt, dx, dy; cx = c.cx, cy = c.cy, kappa = c.kappa, bc = c.bc)
            a, b = b, a
            t += dt
        end
        sqrt(sum(abs2, a - exact_temperature(x, y, tf, id)) / (nx * ny))
    end
    function burgers_error(nx, ny)
        dx, dy = 2 / nx, 1 / ny
        x, y = cell_centers(nx, 2.), cell_centers(ny, 1.)
        u, v = exact_burgers(x, y, 0.)
        un, vn = similar(u), similar(v)
        t = 0.
        while t < tf
            dt = min(N.burgers_stable_timestep(u, v, 0.05, dx, dy), tf - t)
            N.burgers_step!(un, vn, u, v, dt, dx, dy; nu = 0.05)
            u, un = un, u
            v, vn = vn, v
            t += dt
        end
        ue, ve = exact_burgers(x, y, tf)
        sqrt((sum(abs2, u - ue) + sum(abs2, v - ve)) / (nx * ny))
    end
    cases = (
        ("periodic_diffusion", (1.8, 2.2)),
        ("closed_fixed", (1.8, 2.2)),
        ("periodic_combined", (0.8, 1.2)),
        ("burgers", (0.8, 1.2)),
    )
    for (id, (plo, phi)) in cases
        e = [id == "burgers" ? burgers_error(g...) : temperature_error(id, g...) for g in grids]
        p = log2.(e[1:2] ./ e[2:3])
        @test all(>(0), e)             # 誤差が0なら比較が無意味（解析解をそのまま返す誤りを検出）
        @test e[1] > e[2] > e[3]       # 格子を細かくすると単調に減る
        @test all(q -> plo <= q <= phi, p)   # 両区間の観測次数
    end
end
@testset "N07 保存結果と来歴" begin
    @test check_complete()
end
end
