# N05の共通API。周期添字 → 流束 → バッファ検証 → 刻み調整 → 既習課題の入口の順に読みます。
# N01〜N04の更新・境界処理はサブモジュールN01〜N04へ移し、各run.jlの旧入口はここへ委譲します。
# 方針: 抽出前後で数値結果を一致させるため、更新式の演算順序は既習課題のものをそのまま保ちます。

module ThermofluidExercise

# ---------------------------------------------------------------------------
# 共通API
# ---------------------------------------------------------------------------

"""
    check_periodic_index(i, n)

周期添字APIの共通の引数検証（内部用）。`i`,`n`はBoolを除く整数、`n >= 1`、`1 <= i <= n`。
違反は`ArgumentError`。
"""
function check_periodic_index(i, n)
    (i isa Integer && n isa Integer && !(i isa Bool) && !(n isa Bool)) ||
        throw(ArgumentError("添字iと点数nはBoolを除く整数です"))
    n >= 1 || throw(ArgumentError("点数nは1以上です"))
    1 <= i <= n || throw(ArgumentError("添字iは1 <= i <= nです"))
    return nothing
end

"""
    periodic_left_index(i, n)

重複終点を持たない1始まりの周期配列で、左隣の添字を返す。
先頭(`i == 1`)の左隣は末尾`n`、それ以外は`i - 1`。1点(`n == 1`)では`1`。
"""
function periodic_left_index(i, n)
    check_periodic_index(i, n)
    return i == 1 ? n : i - 1
end

"""
    periodic_right_index(i, n)

重複終点を持たない1始まりの周期配列で、右隣の添字を返す。
末尾(`i == n`)の右隣は先頭`1`、それ以外は`i + 1`。1点(`n == 1`)では`1`。
"""
function periodic_right_index(i, n)
    check_periodic_index(i, n)
    return i == n ? 1 : i + 1
end

"""
    finite_scalar(v, name)

有限な実数スカラーか検証する（内部用）。違反は`ArgumentError`。
"""
finite_scalar(v, name) =
    v isa Real && isfinite(v) ? nothing : throw(ArgumentError("$name は有限な実数です"))

"""
    linear_flux(u, speed)

1点の値`u`と速度`speed`から線形流束 f(u) = c u を返す。両方とも有限な実数。
"""
function linear_flux(u, speed)
    finite_scalar(u, "u")
    finite_scalar(speed, "speed")
    return speed * u
end

"""
    burgers_flux(u)

1点の値`u`からBurgers流束 f(u) = u^2/2 を返す。負値も受け付ける有限な実数。
"""
function burgers_flux(u)
    finite_scalar(u, "u")
    return u^2 / 2
end

"""
    validate_buffers(u_new, u_old)

新旧バッファが同じ旧場から安全に更新できるか検証する。
同形状・1始まり・各軸3点以上の浮動小数配列で、記憶領域を共有せず、旧配列の値が有限であること。
新配列の初期値は上書き対象なので検査しない。成功時は`nothing`、違反は書込み前に`ArgumentError`。
"""
function validate_buffers(u_new, u_old)
    (u_new isa AbstractArray{<:AbstractFloat} && u_old isa AbstractArray{<:AbstractFloat}) ||
        throw(ArgumentError("新旧バッファは浮動小数の配列です"))
    Base.require_one_based_indexing(u_new, u_old)
    size(u_new) == size(u_old) || throw(ArgumentError("新旧バッファは同じ形状にしてください"))
    ndims(u_old) >= 1 && all(>=(3), size(u_old)) ||
        throw(ArgumentError("各軸に3点以上が必要です"))
    (u_new === u_old || Base.mightalias(u_new, u_old)) &&
        throw(ArgumentError("新旧バッファが記憶領域を共有しています"))
    all(isfinite, u_old) || throw(ArgumentError("旧配列は有限値にしてください"))
    return nothing
end

"""
    fit_timestep(dt_max, t_final)

刻み上限`dt_max`以下で最終時刻`t_final`へ等間隔に到達する`(; steps, dt)`を返す。
N_t = max(1, ceil(t_f / Δt_max))、Δt = t_f / N_t。`steps`は正整数、`dt`は有限で正の`Float64`。
"""
function fit_timestep(dt_max, t_final)
    for (v, name) in ((dt_max, "dt_max"), (t_final, "t_final"))
        finite_scalar(v, name)
        v isa Bool && throw(ArgumentError("$name にBoolは使えません"))
        v > 0 || throw(ArgumentError("$name は正の値です"))
    end
    count = t_final / dt_max
    isfinite(count) && count < typemax(Int) - 1 ||
        throw(ArgumentError("ステップ数が表現範囲を超えています"))
    steps = max(1, ceil(Int, count))
    dt = Float64(t_final / steps)
    # 割り算の丸めで上限をわずかに超えた場合は、1ステップ増やして上限以下へ戻す。
    if dt > dt_max
        steps += 1
        dt = Float64(t_final / steps)
    end
    isfinite(dt) && dt > 0 || throw(ArgumentError("刻みが表現できません"))
    return (; steps, dt)
end

# ---------------------------------------------------------------------------
# N01: 線形移流（風上・中心差分）と左流入・右流出境界
# ---------------------------------------------------------------------------
module N01

# N01の旧入口はInt配列や同一でない別名(view)も受け付けていたため、共通のvalidate_buffersでは
# 受理範囲が狭まってしまう。契約を保つため、既習N01の検証をこのサブモジュールへ移す。
function validate_step_inputs(u_new, u_old, c, dt, dx)
    u_new === u_old && throw(ArgumentError("新旧で別々のバッファを使ってください"))
    length(u_new) == length(u_old) >= 3 ||
        throw(ArgumentError("二つのバッファを同じ長さの3点以上にしてください"))
    all(isfinite, u_old) || throw(ArgumentError("u_oldの全要素を有限値にしてください"))
    all(isfinite, (c, dt, dx)) || throw(ArgumentError("c、dt、dxを有限値にしてください"))
    c > 0 || throw(ArgumentError("N01では正の移流速度だけを扱います"))
    dt > 0 || throw(ArgumentError("dtは正にしてください"))
    dx > 0 || throw(ArgumentError("dxは正にしてください"))
    return nothing
end

function validate_boundary_inputs(u, left_value)
    length(u) >= 2 || throw(ArgumentError("uには2点以上が必要です"))
    isfinite(left_value) || throw(ArgumentError("left_valueを有限値にしてください"))
    return nothing
end

"""
    upwind_step!(u_new, u_old, c, dt, dx)

正速度の風上差分と陽Euler法で内部点を更新し、`u_new`を返す。両端は旧値を保つ。
"""
function upwind_step!(u_new, u_old, c::Real, dt::Real, dx::Real)
    validate_step_inputs(u_new, u_old, c, dt, dx)
    courant = c * dt / dx
    copyto!(u_new, u_old)
    for i in 2:(length(u_old)-1)
        u_new[i] = u_old[i] - courant * (u_old[i] - u_old[i-1])
    end
    return u_new
end

"""
    centered_step!(u_new, u_old, c, dt, dx)

中心差分と陽Euler法で内部点を更新し、`u_new`を返す。両端は旧値を保つ。
"""
function centered_step!(u_new, u_old, c::Real, dt::Real, dx::Real)
    validate_step_inputs(u_new, u_old, c, dt, dx)
    courant = c * dt / dx
    copyto!(u_new, u_old)
    for i in 2:(length(u_old)-1)
        u_new[i] = u_old[i] - courant / 2 * (u_old[i+1] - u_old[i-1])
    end
    return u_new
end

"""
    apply_boundary!(u; left_value = 1.0)

左端を`left_value`、右端を左隣の値にして入力`u`を返す。2点以上・有限な`left_value`。
"""
function apply_boundary!(u::AbstractVector{<:Real}; left_value::Real=1.0)
    validate_boundary_inputs(u, left_value)
    u[1] = left_value
    u[end] = u[end-1]
    return u
end
end # module N01

# ---------------------------------------------------------------------------
# N02: 保存形Burgers方程式（固定・周期境界）
# ---------------------------------------------------------------------------
module N02
using ..ThermofluidExercise: periodic_left_index, burgers_flux, validate_buffers

function validate_boundary(boundary)
    boundary in (:fixed, :periodic) ||
        throw(ArgumentError("boundaryは:fixedまたは:periodicです"))
    return nothing
end

"""
    nonlinear_upwind_step!(u_new, u_old, dt, dx; boundary = :fixed)

旧配列の流束後退差分と陽Euler法で`u_new`を更新して返す。
固定境界は内部点のみ（両端は旧値）、周期境界は全点を更新する。`u_old`は変更しない。
"""
function nonlinear_upwind_step!(u_new, u_old, dt::Real, dx::Real; boundary=:fixed)
    validate_boundary(boundary)
    (u_new isa AbstractVector && u_old isa AbstractVector) ||
        throw(ArgumentError("新旧バッファは浮動小数の一次元配列です"))
    validate_buffers(u_new, u_old)
    all(v -> isfinite(v) && v > 0, (dt, dx)) ||
        throw(ArgumentError("dtとdxは有限な正の値です"))
    copyto!(u_new, u_old)
    n = length(u_old)
    indices = boundary == :periodic ? (1:n) : (2:(n-1))
    for i in indices
        left = boundary == :periodic ? periodic_left_index(i, n) : i - 1
        u_new[i] = u_old[i] - dt / dx * (burgers_flux(u_old[i]) - burgers_flux(u_old[left]))
    end
    all(isfinite, u_new) || error("更新後に非有限値があります。更新式を確認してください")
    return u_new
end

"""
    apply_boundary!(u; boundary = :fixed)

固定境界では左端を1、右端を左隣の値にする。周期境界では値を変更しない。入力`u`を返す。
"""
function apply_boundary!(u; boundary=:fixed)
    validate_boundary(boundary)
    length(u) >= 3 || throw(ArgumentError("境界配列には3点以上が必要です"))
    if boundary == :fixed
        u[1] = 1.0
        u[end] = u[end-1]
    end
    return u
end
end # module N02

# ---------------------------------------------------------------------------
# N03: 一次元熱拡散（固定温度・断熱）
# ---------------------------------------------------------------------------
module N03
using ..ThermofluidExercise: validate_buffers

positive_finite(v, name) =
    isfinite(v) && v > 0 ? nothing : throw(ArgumentError("$name は有限な正の値です"))

validate_boundary(b) =
    b in (:fixed, :insulated) ? nothing :
    throw(ArgumentError("boundaryは:fixedまたは:insulatedです"))

function check_buffers(u_new, u_old, boundary)
    validate_boundary(boundary)
    (u_new isa AbstractVector && u_old isa AbstractVector) ||
        throw(ArgumentError("新旧バッファは浮動小数の一次元配列です"))
    validate_buffers(u_new, u_old)
    return nothing
end

"""
    apply_boundary!(u_new, u_old, r; boundary = :fixed)

両端だけを更新した`u_new`を返す。固定温度は0、断熱は鏡映ゴースト点で旧温度から更新する。
"""
function apply_boundary!(u_new, u_old, r::Real; boundary=:fixed)
    check_buffers(u_new, u_old, boundary)
    positive_finite(r, "r")
    if boundary == :fixed
        u_new[1] = 0.0
        u_new[end] = 0.0
    else
        u_new[1] = u_old[1] + 2r * (u_old[2] - u_old[1])
        u_new[end] = u_old[end] + 2r * (u_old[end-1] - u_old[end])
    end
    return u_new
end

"""
    diffusion_step!(u_new, u_old, dt, dx, diffusivity; boundary = :fixed)

中心差分と陽Euler法で内部点を更新し、続けて両端の境界を適用した`u_new`を返す。
安定上限を超える刻みも実験用に受け付ける。
"""
function diffusion_step!(
    u_new, u_old, dt::Real, dx::Real, diffusivity::Real; boundary=:fixed,
)
    check_buffers(u_new, u_old, boundary)
    for (value, name) in ((dt, "dt"), (dx, "dx"), (diffusivity, "diffusivity"))
        positive_finite(value, name)
    end
    r = diffusivity * dt / dx^2
    positive_finite(r, "実効Fourier数")
    for i in 2:(length(u_old)-1)
        u_new[i] = u_old[i] + r * (u_old[i+1] - 2u_old[i] + u_old[i-1])
    end
    apply_boundary!(u_new, u_old, r; boundary)
    return u_new
end
end # module N03

# ---------------------------------------------------------------------------
# N04: 周期一次元移流拡散（モデルは呼出し側が引数で渡す）
# ---------------------------------------------------------------------------
module N04
using ..ThermofluidExercise:
    periodic_left_index, periodic_right_index, linear_flux, burgers_flux, validate_buffers

finite_real(v, name) =
    v isa Real && isfinite(v) ? nothing : throw(ArgumentError("$name は有限な実数です"))

function nonnegative_finite(v, name)
    finite_real(v, name)
    v >= 0 || throw(ArgumentError("$name は非負です"))
end

function positive_finite(v, name)
    finite_real(v, name)
    v > 0 || throw(ArgumentError("$name は正です"))
end

validate_model(m) =
    m in (:linear, :nonlinear) ? nothing :
    throw(ArgumentError("modelは:linearまたは:nonlinearです（model=$m）"))

function validate_speed(model, speed)
    finite_real(speed, "speed")
    model == :linear && positive_finite(speed, "線形speed")
end

function validate_safety(safety)
    positive_finite(safety, "safety")
    safety <= 1 || throw(ArgumentError("safetyは0より大きく1以下です"))
end

function maximum_speed(u, model, advection, speed)
    !advection && return 0.0
    if model == :nonlinear
        minimum(u) >= 0 ||
            throw(ArgumentError("非線形移流は非負速度の範囲だけに対応します"))
        return maximum(u)
    end
    return speed
end

function stability_numbers(a, dt, dx, d)
    cfl = a * (dt / dx)
    fo = d * (dt / dx) / dx
    all(isfinite, (cfl, fo)) || throw(ArgumentError("CFL・Foが有限ではありません"))
    cfl + 2fo <= 1 + 32eps(Float64) ||
        throw(ArgumentError("合成安定条件CFL+2Fo<=1を超えています: CFL=$cfl, Fo=$fo"))
    return cfl, fo
end

"""
    advective_flux(u; model, speed = 1.0)

選択モデルの保存形流束。線形は共通の`linear_flux`、非線形（非負のu）は共通の`burgers_flux`。
"""
function advective_flux(u; model, speed=1.0)
    validate_model(model)
    finite_real(u, "u")
    validate_speed(model, speed)
    if model == :linear
        return linear_flux(u, speed)
    else
        nonnegative_finite(u, "非線形移流のu")
        return burgers_flux(u)
    end
end

"""
    stable_timestep(max_speed, dx, diffusivity; safety = 0.8)

合成安定条件 CFL + 2Fo <= 1 から求めた上限に`safety`を掛けた刻みを返す。
"""
function stable_timestep(max_speed, dx, diffusivity; safety=0.8)
    nonnegative_finite(max_speed, "max_speed")
    positive_finite(dx, "dx")
    nonnegative_finite(diffusivity, "diffusivity")
    validate_safety(safety)
    max_speed > 0 || diffusivity > 0 ||
        throw(ArgumentError("移流速度とdiffusivityがともに0では刻みを定義できません"))
    return safety / (max_speed / dx + 2diffusivity / dx^2)
end

"""
    advection_diffusion_step!(u_new, u_old, dt, dx, diffusivity; model, advection = true, speed = 1.0)

周期の左右隣接を使い、移流と拡散を同じ旧配列から評価して全点を更新した`u_new`を返す。
`model`は既定値を持たず、呼出し側（課題のSELECTED_MODEL）から必ず渡す。
"""
function advection_diffusion_step!(
    u_new, u_old, dt, dx, diffusivity; model, advection=true, speed=1.0,
)
    validate_model(model)
    validate_speed(model, speed)
    advection isa Bool || throw(ArgumentError("advectionはBoolです"))
    (u_new isa AbstractVector && u_old isa AbstractVector) ||
        throw(ArgumentError("新旧配列は浮動小数の一次元配列です"))
    validate_buffers(u_new, u_old)
    positive_finite(dt, "dt")
    positive_finite(dx, "dx")
    nonnegative_finite(diffusivity, "diffusivity")
    stability_numbers(maximum_speed(u_old, model, advection, speed), dt, dx, diffusivity)

    n = length(u_old)
    for i in eachindex(u_old)
        im = periodic_left_index(i, n)
        ip = periodic_right_index(i, n)
        adv = advection ?
              dt / dx * (advective_flux(u_old[i]; model, speed) -
            advective_flux(u_old[im]; model, speed)) : 0.0
        dif = diffusivity * dt / dx^2 * (u_old[ip] - 2u_old[i] + u_old[im])
        u_new[i] = u_old[i] - adv + dif
    end
    return u_new
end
end # module N04

include("N06Advection.jl")
include("N07Transport.jl")
include("N08N09Elliptic.jl")
end
