# N02: 流束、周期左隣添字、流束差分のTODO 3か所を実装する。
# 読む順: 流束 → 添字 → 1ステップ → 境界 → simulate → main。
module N02NonlinearAdvection
include("provided_support.jl")
export burgers_flux, periodic_left_index, nonlinear_upwind_step!, apply_boundary!, simulate, main

"""有限なスカラーのBurgers流束。流束自体は負の値にも定義される。"""
function burgers_flux(u::Real)
    isfinite(u) || throw(ArgumentError("流束の入力を有限値にしてください"))
    # TODO(N02): 保存形Burgers方程式の流束を返す。
    error("未実装 N02: burgers_fluxを実装してください")
end

"""1始まりの周期配列で左隣の添字を返す。"""
function periodic_left_index(i::Integer, n::Integer)
    (i isa Bool || n isa Bool || n < 1 || !(1 <= i <= n)) &&
        throw(ArgumentError("添字は1 <= i <= n、点数はn >= 1の整数です（Bool不可）"))
    # TODO(N02): 先頭の左隣が末尾になるようにする。
    error("未実装 N02: periodic_left_indexを実装してください")
end

"""旧配列だけを読み、流束後退差分と陽Eulerで新配列へ書く。
固定境界では端点を保持し、周期境界では全点を更新する。
新旧は独立した同長3点以上の浮動小数配列。境界適用は呼出し側で行う。
"""
function nonlinear_upwind_step!(u_new, u_old, dt::Real, dx::Real; boundary=:fixed)
    validate_step_inputs(u_new, u_old, dt, dx, boundary)
    copyto!(u_new, u_old)
    indices = boundary == :periodic ? (1:length(u_old)) : (2:length(u_old)-1)
    for i in indices
        left = boundary == :periodic ? periodic_left_index(i, length(u_old)) : i - 1
        # TODO(N02): 古い値の流束の差で更新する。u_oldは変更しない。
        error("未実装 N02: nonlinear_upwind_step!の更新式を実装してください")
    end
    all(isfinite, u_new) || error("更新後に非有限値があります。更新式を確認してください")
    return u_new
end

"""固定は左1・右隣接点コピー、周期は配列を変更しない。"""
function apply_boundary!(u; boundary=:fixed)
    validate_boundary(boundary)
    length(u) >= 3 || throw(ArgumentError("境界配列には3点以上が必要です"))
    if boundary == :fixed
        u[1] = 1.0
        u[end] = u[end-1]
    end
    return u
end

"""独立したx,u0,uと格子・時間刻み・初期/最終診断を返す。
周期公式条件はnx=80を明示する。初期最大速度から刻みを定め、最終時刻に合わせる。
CFL超過でも計算を続ける。発散して非有限値になった場合は停止する。
"""
function simulate(; boundary=:fixed, nx::Integer=81, cfl::Real=0.5, t_final::Real=1.0)
    validate_simulation_inputs(boundary, nx, cfl, t_final)
    dx = 2.0 / (boundary == :fixed ? nx-1 : nx)
    x = [j * dx for j in 0:nx-1]
    u0 = [0.5 <= xi <= 1.0 ? 2.0 : 1.0 for xi in x]
    nominal_dt = cfl * dx / maximum(abs, u0)
    steps = ceil(Int, t_final / nominal_dt)
    dt = t_final / steps
    u_old = copy(u0)
    u_new = similar(u_old)
    max_cfl = maximum(abs, u_old) * dt / dx
    for _ in 1:steps
        nonlinear_upwind_step!(u_new, u_old, dt, dx; boundary)
        apply_boundary!(u_new; boundary)
        all(isfinite, u_new) || error("更新後に非有限値があります")
        max_cfl = max(max_cfl, maximum(abs, u_new) * dt / dx)
        u_old, u_new = u_new, u_old
    end
    initial_minimum, initial_maximum = extrema(u0)
    low, high = extrema(u_old)
    return (; x, u0, u=u_old, dx, dt, steps, max_cfl, t_final,
        initial_minimum, initial_maximum, minimum=low, maximum=high,
        overshoot=max(high-initial_maximum,0.0), undershoot=max(initial_minimum-low,0.0),
        initial_sum=sum(u0), final_sum=sum(u_old), sum_change=sum(u_old)-sum(u0))
end

"""固定81点・周期80点の公式計算と3出力を生成する。"""
function main(; output_dir::AbstractString=DEFAULT_OUTPUT_DIR, cfl::Real=0.5, t_final::Real=1.0)
    fixed = simulate(; boundary=:fixed, nx=81, cfl, t_final)
    periodic = simulate(; boundary=:periodic, nx=80, cfl, t_final)
    write_outputs(output_dir, fixed, periodic, cfl, t_final)
    println("N02の出力を書き込みました: $(abspath(output_dir))")
    return (; fixed, periodic)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
end
