# N04: SELECTED_MODELを選び、選択流束・合成刻み・周期更新の3か所を実装する。
# 読む順: 流束 → 刻み → 1ステップ → simulate → main。入力検証・解析解・描画は提供済み。
module N04AdvectionDiffusion
const SELECTED_MODEL = :unselected
include("provided_support.jl")
export SELECTED_MODEL, advective_flux, stable_timestep, advection_diffusion_step!,
    initial_condition, analytic_solution, conserved_integral, simulate, main

"""選択モデルの保存形移流流束。"""
function advective_flux(u; model=SELECTED_MODEL, speed=1.0)
    validate_model(model); finite_real(u,"u"); validate_speed(model,speed)
    if model == :linear
        # TODO(N04): 線形を選んだ場合だけ実装する。
        error("未実装 N04: 線形流束を実装してください")
    else
        nonnegative_finite(u,"非線形移流のu")
        # TODO(N04): 非線形を選んだ場合だけ実装する。
        error("未実装 N04: 非線形流束を実装してください")
    end
end

"""C+2Foの合成安定条件を満たす刻み。"""
function stable_timestep(max_speed,dx,diffusivity; safety=0.8)
    validate_timestep_inputs(max_speed,dx,diffusivity,safety)
    # TODO(N04): 合成条件から刻みを返す。
    error("未実装 N04: stable_timestepを実装してください")
end

"""u_oldから周期全点の次の値を計算し、u_newへ書き込む。"""
function advection_diffusion_step!(u_new,u_old,dt,dx,diffusivity;
    model=SELECTED_MODEL,advection=true,speed=1.0)
    validate_step(u_new,u_old,dt,dx,diffusivity,model,advection,speed)
    for i in eachindex(u_old)
        # TODO(N04): 周期の左右隣接を使い、移流と拡散を同じ旧配列から足す。
        error("未実装 N04: 周期更新を実装してください")
    end
    return u_new
end

"""指定最終時刻まで計算し、初期・最終配列、積分履歴、全時刻の診断を返す。
dtを指定した場合は、その値を上限に最終時刻へ合わせて刻みを調整する。
"""
function simulate(;model=SELECTED_MODEL,advection=true,nx=80,diffusivity=0.1,
    speed=1.0,safety=0.8,t_final=1.0,initial=:pulse,dt=nothing)
    validate_simulation(model,advection,nx,diffusivity,speed,safety,t_final,initial)
    dx=2.0/nx
    x=[j*dx for j in 0:nx-1]
    u0=initial_condition(x;model,initial,diffusivity,speed)
    a=maximum_speed(u0,model,advection,speed)
    nominal_dt=stable_timestep(a,dx,diffusivity;safety)
    positive_finite(nominal_dt,"nominal_dt")
    if !isnothing(dt)
        positive_finite(dt,"dt")
        stability_numbers(a,dt,dx,diffusivity)
        nominal_dt=dt
    end
    count=t_final/nominal_dt
    isfinite(count) && count < typemax(Int)-1 || throw(ArgumentError("ステップ数が表現範囲を超えています"))
    steps=max(1,ceil(Int,count))
    dt=t_final/steps
    positive_finite(dt,"実効dt")
    cfl,fo=stability_numbers(a,dt,dx,diffusivity)
    u_old=copy(u0); u_new=similar(u_old)
    times=collect(range(0.,t_final;length=steps+1))
    integral_history=Vector{Float64}(undef,steps+1)
    initial_integral=conserved_integral(u0,dx)
    integral_history[1]=initial_integral
    low,high=extrema(u0)
    max_cfl=cfl
    for step in 1:steps
        a=maximum_speed(u_old,model,advection,speed)
        cfl,fo=stability_numbers(a,dt,dx,diffusivity)
        advection_diffusion_step!(u_new,u_old,dt,dx,diffusivity;model,advection,speed)
        ensure_finite(u_new,model,step,times[step+1],cfl,fo)
        integral=dx*sum(u_new)
        ensure_finite((integral,),model,step,times[step+1],cfl,fo)
        integral_history[step+1]=integral
        low=min(low,minimum(u_new)); high=max(high,maximum(u_new))
        a=maximum_speed(u_new,model,advection,speed)
        cfl,fo=stability_numbers(a,dt,dx,diffusivity)
        max_cfl=max(max_cfl,cfl)
        u_old,u_new=u_new,u_old
    end
    final_integral=integral_history[end]
    integral_change=final_integral-initial_integral
    max_stability_number=max_cfl+2fo
    ensure_finite((integral_change,low,high,max_stability_number),model,steps,t_final,max_cfl,fo)
    return (;x,u0,u=u_old,model,advection,initial,speed,diffusivity,dx,dt,steps,t_final,
        requested_safety=safety,times,integral_history,initial_integral,final_integral,integral_change,
        minimum=low,maximum=high,max_cfl,fo,max_stability_number)
end

"""選択モデルの3条件と解析モードの格子収束をresults/<model>/へ出力する。"""
function main(;model=SELECTED_MODEL,output_dir=DEFAULT_OUTPUT_DIR)
    validate_model(model)
    combined=simulate(;model)
    advection_only=simulate(;model,diffusivity=0.,dt=combined.dt)
    diffusion_only=simulate(;model,advection=false,dt=combined.dt)
    convergence=convergence_results(model)
    write_outputs(output_dir,advection_only,diffusion_only,combined,convergence)
    println("N04の出力を書き込みました: $(abspath(joinpath(output_dir,string(model))))")
    return (;advection_only,diffusion_only,combined,convergence)
end

if abspath(PROGRAM_FILE)==@__FILE__
    main()
end
end
