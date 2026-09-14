# N01: 初期条件、風上差分、中心差分のTODO 3か所を実装する。
# 読む順: 初期条件 → 1ステップ更新 → 境界条件 → simulate → 診断 → main。
# 自分の確認はtests.jl、記録はlearning_log.md。実行結果はresults/へ保存する。
# 入力検証・描画・保存の詳細はprovided_support.jlを参照する。

module N01LinearAdvection

include("provided_support.jl")

export rectangular_initial_condition, apply_boundary!, upwind_step!, centered_step!
export simulate, write_summary, make_plots, main

# === 学生が実装する3つの関数 ===

"""
    rectangular_initial_condition(x; base, plateau, plateau_start, plateau_end)

N01の移流実験で使う矩形状の初期分布を作る。
各座標を読み、矩形領域の内側を`plateau`、外側を`base`にする。

# 引数

- `x`: 狭義単調増加する座標の配列。
- `base`: 矩形領域の外側の値。
- `plateau`: 矩形領域の内側の値。
- `plateau_start`: 矩形領域を始める座標。
- `plateau_end`: 矩形領域を終える座標。

# 戻り値

`x`と同じ長さを持つ初期値の配列を返す。`x`は変更しない。
"""
function rectangular_initial_condition(
    x::AbstractVector{<:Real};
    base::Real = 1.0,
    plateau::Real = 2.0,
    plateau_start::Real = 0.5,
    plateau_end::Real = 1.0,
)
    validate_initial_condition_inputs(x, base, plateau, plateau_start, plateau_end)

    #=
    TODO(N01):
    矩形状の初期分布を実装する。
    各座標xiがplateau_start <= xi <= plateau_endを満たすか判定する。
    =#
    return fill(float(base), length(x))
end

"""
    upwind_step!(u_new, u_old, c, dt, dx)

正の移流速度に対する風上差分と陽Euler法で1ステップ進める。
`u_old`だけを読み、計算結果を`u_new`へ書き換える。`u_old`は変更しない。

# 引数

- `u_new`: 新しい時刻の値を書き込む配列。
- `u_old`: 現在時刻の値を読む配列。
- `c`: 正の移流速度。
- `dt`: 時間刻み。
- `dx`: 格子間隔。

# 戻り値

書き換えた`u_new`を返す。
"""
function upwind_step!(u_new, u_old, c::Real, dt::Real, dx::Real)
    validate_step_inputs(u_new, u_old, c, dt, dx)
    # courantは1ステップで進む格子幅の割合を表すCFL数です。
    courant = c * dt / dx
    # copyto!で古い値を複製し、内部を更新する間も両端点の値を保ちます。
    copyto!(u_new, u_old)
    for i in 2:(length(u_old) - 1)
        #=
        TODO(N01):
        風上差分と陽Euler法による1ステップの更新を実装する。
        授業ページの式を、iとi - 1の添字を使ってコードへ写す。
        =#
        u_new[i] = u_old[i]
    end
    return u_new
end

"""
    centered_step!(u_new, u_old, c, dt, dx)

意図的に不安定な中心差分と陽Euler法で1ステップ進める。
`u_old`だけを読み、計算結果を`u_new`へ書き換える。`u_old`は変更しない。

# 引数

- `u_new`: 新しい時刻の値を書き込む配列。
- `u_old`: 現在時刻の値を読む配列。
- `c`: 正の移流速度。
- `dt`: 時間刻み。
- `dx`: 格子間隔。

# 戻り値

書き換えた`u_new`を返す。
"""
function centered_step!(u_new, u_old, c::Real, dt::Real, dx::Real)
    validate_step_inputs(u_new, u_old, c, dt, dx)
    # courantは1ステップで進む格子幅の割合を表すCFL数です。
    courant = c * dt / dx
    # copyto!で古い値を複製し、内部を更新する間も両端点の値を保ちます。
    copyto!(u_new, u_old)
    for i in 2:(length(u_old) - 1)
        #=
        TODO(N01):
        中心差分と陽Euler法による1ステップの更新を実装する。
        授業ページの式を、i - 1とi + 1の添字を使ってコードへ写す。
        =#
        u_new[i] = u_old[i]
    end
    return u_new
end

# === 境界条件と時間発展の流れ ===

"""
    apply_boundary!(u; left_value)

1ステップ更新後の配列`u`へ境界条件を適用する。
左端を`left_value`に固定し、右端を左隣と同じ値にしてゼロ勾配を表す。

# 引数

- `u`: 境界値を書き換える配列。
- `left_value`: 左端へ設定する値。

# 戻り値

境界値を書き換えた`u`を返す。
"""
function apply_boundary!(u::AbstractVector{<:Real}; left_value::Real = 1.0)
    validate_boundary_inputs(u, left_value)
    u[1] = left_value
    u[end] = u[end - 1]
    return u
end

"""
    simulate(; scheme, nx, c, cfl, t_final)

指定した差分法でN01の時間発展を計算する。
初期条件を作り、選択した1ステップ関数と境界条件を各ステップで順に適用する。

# 引数

- `scheme`: `:upwind`または`:centered`。
- `nx`: 0から2までに置く格子点数。
- `c`: 正の移流速度。
- `cfl`: 1ステップで進む格子幅の割合。
- `t_final`: 計算する最終時刻。

# 戻り値

次のフィールドを持つ名前付きタプルを返す。

- `x`: 座標。
- `u0`: 初期値。
- `u`: 最終時刻の値。
- `dx`, `dt`, `steps`, `cfl`: 格子幅、時間刻み、ステップ数、実効CFL。
- `minimum`, `maximum`: 最終値の最小値と最大値。
"""
function simulate(;
    scheme,
    nx::Integer = 81,
    c::Real = 1.0,
    cfl::Real = 0.5,
    t_final::Real = 0.5,
)
    validate_simulation_inputs(scheme, nx, c, cfl, t_final)

    x = collect(range(0.0, 2.0; length = nx))
    dx = x[2] - x[1]
    nominal_dt = cfl * dx / c
    steps = ceil(Int, t_final / nominal_dt)
    dt = t_final / steps
    actual_cfl = c * dt / dx

    u0 = rectangular_initial_condition(x)
    u_old = copy(u0)
    u_new = similar(u_old)
    # `scheme`に応じて、時間ループで呼ぶ1ステップ関数を選びます。
    step! = scheme === :upwind ? upwind_step! : centered_step!

    # 各ステップで新しい値を計算し、境界条件を適用してから二つのバッファを交換する。
    for _ in 1:steps
        step!(u_new, u_old, c, dt, dx)
        apply_boundary!(u_new)
        u_old, u_new = u_new, u_old
    end

    return (
        x = x,
        u0 = u0,
        u = u_old,
        dx = dx,
        dt = dt,
        steps = steps,
        cfl = actual_cfl,
        minimum = minimum(u_old),
        maximum = maximum(u_old),
    )
end

"""一つの差分法についてTOMLへ書き出す診断量を作る。"""
function summary_section(scheme::String, result)
    initial_minimum, initial_maximum = extrema(result.u0)
    overshoot = max(result.maximum - initial_maximum, 0.0)
    undershoot = max(initial_minimum - result.minimum, 0.0)
    tolerance = 100eps(Float64) * max(abs(initial_minimum), abs(initial_maximum), 1.0)
    return Dict(
        "scheme" => scheme,
        "cfl" => result.cfl,
        "dt" => result.dt,
        "steps" => result.steps,
        "minimum" => result.minimum,
        "maximum" => result.maximum,
        "overshoot" => overshoot,
        "undershoot" => undershoot,
        "overshoot_occurred" => overshoot > tolerance,
        "undershoot_occurred" => undershoot > tolerance,
    )
end

"""N01の二つの比較計算を実行し、公式出力をすべて書き出す。"""
function main(;
    output_dir::AbstractString = DEFAULT_OUTPUT_DIR,
    nx::Integer = 81,
    c::Real = 1.0,
    cfl::Real = 0.5,
    t_final::Real = 0.5,
)
    upwind = simulate(; scheme = :upwind, nx, c, cfl, t_final)
    centered = simulate(; scheme = :centered, nx, c, cfl, t_final)
    summary_path = write_summary(output_dir, upwind, centered)
    plot_paths = make_plots(output_dir, upwind, centered)
    println("N01の出力を書き込みました: $(abspath(output_dir))")
    return (; upwind, centered, summary_path, plot_paths)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end
