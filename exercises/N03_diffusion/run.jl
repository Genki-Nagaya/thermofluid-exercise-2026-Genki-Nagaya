# N03: 内部点更新、断熱端更新、熱量積分のTODO 3か所を実装する。
# 読む順: 初期条件 → 境界 → 1ステップ → 熱量 → simulate → main。
module N03Diffusion
include("provided_support.jl")
export initial_condition, analytic_solution, apply_boundary!, diffusion_step!, thermal_content,
    simulate, main, stability_experiment

"""無次元領域[0,2]上の矩形、または境界に適合する解析モード。"""
function initial_condition(x; initial=:pulse, boundary=:fixed)
    validate_boundary(boundary)
    validate_initial(initial)
    all(isfinite,x) || throw(ArgumentError("座標を有限値にしてください"))
    initial == :mode && return analytic_solution(x,0.;boundary)
    return [0.5 <= xi <= 1.0 ? 1.0 : 0.0 for xi in x]
end

"""解析モード専用の解。矩形初期値の解析解ではない。"""
function analytic_solution(x,t::Real; boundary=:fixed,diffusivity::Real=0.1)
    validate_boundary(boundary)
    positive_finite(diffusivity,"diffusivity")
    isfinite(t) && t >= 0 || throw(ArgumentError("tは有限な非負値です"))
    all(isfinite,x) || throw(ArgumentError("座標を有限値にしてください"))
    amplitude=exp(-diffusivity*(pi/2)^2*t)
    return boundary == :fixed ? amplitude .* sin.(pi .* x ./ 2) :
        0.5 .+ 0.5 .* amplitude .* cos.(pi .* x ./ 2)
end

"""固定端は0。断熱端は旧配列の鏡映ゴースト点から更新する。"""
function apply_boundary!(u_new,u_old,r::Real;boundary=:fixed)
    validate_buffers(u_new,u_old,boundary)
    positive_finite(r,"r")
    if boundary == :fixed
        u_new[1]=0.0
        u_new[end]=0.0
    else
        # TODO(N03): 旧配列を用い、半セルに対応する係数で両端を更新する。
        error("未実装 N03: 断熱端の更新を実装してください")
    end
    return u_new
end

"""空間中心差分・陽Euler。境界もここで更新し、旧配列は変更しない。"""
function diffusion_step!(u_new,u_old,dt::Real,dx::Real,diffusivity::Real;boundary=:fixed)
    validate_buffers(u_new,u_old,boundary)
    for (value,name) in ((dt,"dt"),(dx,"dx"),(diffusivity,"diffusivity"))
        positive_finite(value,name)
    end
    r=diffusivity*dt/dx^2
    positive_finite(r,"実効Fourier数")
    for i in 2:length(u_old)-1
        # TODO(N03): 旧配列だけから内部点の次時刻値を計算する。
        error("未実装 N03: 内部点の更新を実装してください")
    end
    apply_boundary!(u_new,u_old,r;boundary)
    return u_new
end

"""端点に半分の重みを付けた無次元熱量H。"""
function thermal_content(u,dx::Real)
    u isa AbstractVector{<:Real} && length(u)>=3 || throw(ArgumentError("熱量配列は3点以上です"))
    Base.require_one_based_indexing(u)
    all(isfinite,u) || throw(ArgumentError("温度を有限値にしてください"))
    positive_finite(dx,"dx")
    # TODO(N03): 台形則の積分を返す。
    error("未実装 N03: thermal_contentを実装してください")
end

"""指定時刻へ刻みを合わせる。foは実効値、requested_foは指定値。
温度全履歴は保持せず、熱量履歴と初期・最終温度、全時刻の極値を返す。
安定条件超過や負値も許容し、非有限値はstep・t・foを示して停止する。
"""
function simulate(;boundary=:fixed,nx::Integer=81,diffusivity::Real=0.1,
    fo::Real=0.4,t_final::Real=1.0,initial=:pulse)
    validate_simulation_inputs(boundary,nx,diffusivity,fo,t_final,initial)
    dx=2.0/(nx-1)
    x=[j*dx for j in 0:nx-1]
    u0=initial_condition(x;initial,boundary)
    # 正弦の浮動小数丸めも含め、固定端を正確に0へ揃える。
    if boundary == :fixed
        u0[1]=0.; u0[end]=0.
    end
    nominal_dt=fo*dx^2/diffusivity
    positive_finite(nominal_dt,"nominal_dt")
    count=t_final/nominal_dt
    isfinite(count) && count < typemax(Int) || throw(ArgumentError("ステップ数が表現範囲を超えています"))
    steps=max(1,ceil(Int,count))
    dt=t_final/steps
    effective_fo=diffusivity*dt/dx^2
    positive_finite(effective_fo,"実効Fourier数")
    u_old=copy(u0); u_new=similar(u_old)
    initial_heat=thermal_content(u0,dx)
    initial_minimum,initial_maximum=extrema(u0)
    low,high=initial_minimum,initial_maximum
    ensure_finite((initial_heat,low,high),0,0.,effective_fo)
    times=collect(range(0.,t_final;length=steps+1))
    heat_history=Vector{Float64}(undef,steps+1)
    heat_history[1]=initial_heat
    for step in 1:steps
        diffusion_step!(u_new,u_old,dt,dx,diffusivity;boundary)
        ensure_finite(u_new,step,times[step+1],effective_fo)
        heat=thermal_content(u_new,dx)
        ensure_finite((heat,),step,times[step+1],effective_fo)
        heat_history[step+1]=heat
        low=min(low,minimum(u_new)); high=max(high,maximum(u_new))
        u_old,u_new=u_new,u_old
    end
    final_heat=heat_history[end]
    heat_change=final_heat-initial_heat
    ensure_finite((heat_change,),steps,t_final,effective_fo)
    return (;x,u0,u=u_old,dx,dt,steps,fo=effective_fo,requested_fo=fo,t_final,
        diffusivity,times,heat_history,initial_heat,final_heat,heat_change,
        initial_minimum,initial_maximum,minimum=low,maximum=high)
end

"""両境界の標準計算と解析モードの格子収束から公式4出力を生成する。"""
function main(;output_dir::AbstractString=DEFAULT_OUTPUT_DIR)
    fixed=simulate()
    insulated=simulate(;boundary=:insulated)
    convergence=convergence_results()
    write_outputs(output_dir,fixed,insulated,convergence)
    println("N03の出力を書き込みました: $(abspath(output_dir))")
    return (;fixed,insulated,convergence)
end

"""任意実験。固定温度で安定・不安定条件を同時刻まで比較する。"""
function stability_experiment(;output_dir::AbstractString,fo::Real=0.6,t_final::Real=0.1)
    isfinite(fo) && fo>0.5 || throw(ArgumentError("比較用foは有限で0.5より大きくしてください"))
    stable=simulate(;fo=0.4,t_final)
    unstable=simulate(;fo,t_final)
    unstable.fo>0.5 || throw(ArgumentError("最終時刻への調整後の実効foが0.5以下です。安定条件超過の比較になるよう条件を調整してください"))
    write_stability_outputs(output_dir,stable,unstable)
    return (;stable,unstable)
end

if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
end
