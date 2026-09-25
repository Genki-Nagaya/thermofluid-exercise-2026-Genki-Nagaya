# 入力検証・解析解・診断・描画・保存の提供部分。学生の編集対象ではありません。
using Plots, TOML
const DEFAULT_OUTPUT_DIR=joinpath(@__DIR__,"results")
finite_real(v,name) = v isa Real && isfinite(v) ? nothing : throw(ArgumentError("$name は有限な実数です"))
function nonnegative_finite(v,name)
    finite_real(v,name)
    v>=0 || throw(ArgumentError("$name は非負です"))
end
function positive_finite(v,name)
    finite_real(v,name)
    v>0 || throw(ArgumentError("$name は正です"))
end
validate_model(m) = m in (:linear,:nonlinear) ? nothing : throw(ArgumentError("SELECTED_MODELを:linearまたは:nonlinearに設定してください（model=$m）"))
function validate_speed(model,speed)
    finite_real(speed,"speed")
    model==:linear && positive_finite(speed,"線形speed")
end
function validate_safety(safety)
    positive_finite(safety,"safety")
    safety<=1 || throw(ArgumentError("safetyは0より大きく1以下です"))
end
function validate_timestep_inputs(a,dx,d,safety)
    nonnegative_finite(a,"max_speed"); positive_finite(dx,"dx")
    nonnegative_finite(d,"diffusivity"); validate_safety(safety)
    a>0 || d>0 || throw(ArgumentError("移流速度とdiffusivityがともに0では刻みを定義できません"))
end
function checked_timestep(dt)
    positive_finite(dt,"合成条件のdt（係数・dxの表現範囲を確認）")
    return dt
end
function maximum_speed(u,model,advection,speed)
    !advection && return 0.
    if model==:nonlinear
        minimum(u)>=0 || throw(ArgumentError("非線形移流は非負速度の範囲だけに対応します"))
        return maximum(u)
    end
    return speed
end
function stability_numbers(a,dt,dx,d)
    cfl=a*(dt/dx); fo=d*(dt/dx)/dx
    all(isfinite,(cfl,fo)) || throw(ArgumentError("CFL・Foが有限ではありません"))
    cfl+2fo<=1+32eps(Float64) || throw(ArgumentError("合成安定条件CFL+2Fo<=1を超えています: CFL=$cfl, Fo=$fo"))
    return cfl,fo
end
function validate_step(a,b,dt,dx,d,model,advection,speed)
    validate_model(model); validate_speed(model,speed)
    advection isa Bool || throw(ArgumentError("advectionはBoolです"))
    all(u->u isa AbstractVector{<:AbstractFloat},(a,b)) || throw(ArgumentError("新旧配列は浮動小数の一次元配列です"))
    Base.require_one_based_indexing(a,b)
    length(a)==length(b)>=3 || throw(ArgumentError("新旧配列は同長で3点以上です"))
    Base.mightalias(a,b) && throw(ArgumentError("新旧配列が重なっています"))
    all(isfinite,b) || throw(ArgumentError("旧配列は有限値にしてください"))
    positive_finite(dt,"dt"); positive_finite(dx,"dx"); nonnegative_finite(d,"diffusivity")
    stability_numbers(maximum_speed(b,model,advection,speed),dt,dx,d)
end
function validate_simulation(model,advection,nx,d,speed,safety,t,initial)
    validate_model(model); validate_speed(model,speed); validate_safety(safety)
    advection isa Bool || throw(ArgumentError("advectionはBoolです"))
    nx isa Integer && !(nx isa Bool) && nx>=3 || throw(ArgumentError("nxは3以上の整数です（Bool不可）"))
    nonnegative_finite(d,"diffusivity"); positive_finite(t,"t_final")
    initial in (:pulse,:mode) || throw(ArgumentError("initialは:pulseまたは:modeです"))
end
function ensure_finite(values,model,step,t,cfl,fo)
    all(isfinite,values) || error("非有限値のため停止: model=$model, step=$step, t=$t, CFL=$cfl, Fo=$fo")
end
"""周期[0,2)の矩形、または選択モデルの解析モードの初期値。"""
function initial_condition(x;model=SELECTED_MODEL,initial=:pulse,diffusivity=0.1,speed=1.0)
    validate_model(model); validate_speed(model,speed); nonnegative_finite(diffusivity,"diffusivity")
    all(isfinite,x) || throw(ArgumentError("xは有限値です"))
    initial in (:pulse,:mode) || throw(ArgumentError("initialは:pulseまたは:modeです"))
    initial==:mode && return analytic_solution(x,0.;model,diffusivity,speed)
    return [0.5<=xi<=1. ? 2. : 1. for xi in x]
end
"""滑らかな解析モード専用。矩形の解析解ではない。非線形の平均速度は1。"""
function analytic_solution(x,t;model=SELECTED_MODEL,diffusivity=0.1,speed=1.0)
    validate_model(model); validate_speed(model,speed)
    nonnegative_finite(diffusivity,"diffusivity"); nonnegative_finite(t,"t")
    all(isfinite,x) || throw(ArgumentError("xは有限値です"))
    if model==:linear
        return 1 .+ 0.5exp(-diffusivity*pi^2*t).*sin.(pi.*(x.-speed*t))
    end
    positive_finite(diffusivity,"非線形解析モードのdiffusivity")
    # max |u-1| at t=0 is 2D*k*A/sqrt(1-A^2), A=0.5.
    diffusivity*pi/sqrt(0.75)<=1 || throw(ArgumentError("非線形解析モードで負速度を生むdiffusivityです"))
    a=0.5exp(-diffusivity*pi^2*t); theta=pi.*(x.-t)
    return 1 .+ 2diffusivity*pi*a.*sin.(theta)./(1 .+ a.*cos.(theta))
end
"""周期等幅セルの離散積分。終点を重複させず、半重みを付けない。"""
function conserved_integral(u,dx)
    u isa AbstractVector{<:Real} && length(u)>=3 || throw(ArgumentError("積分配列は3点以上の実数ベクトルです"))
    Base.require_one_based_indexing(u)
    all(isfinite,u) || throw(ArgumentError("積分配列は有限値です"))
    positive_finite(dx,"dx")
    result=dx*sum(u)
    finite_real(result,"離散積分")
    return result
end
function convergence_results(model)
    nx=[40,80,160]
    runs=[simulate(;model,nx=n,initial=:mode) for n in nx]
    errors=[maximum(abs.(r.u-analytic_solution(r.x,r.t_final;model))) for r in runs]
    for (r,error) in zip(runs,errors)
        ensure_finite((error,),model,r.steps,r.t_final,r.max_cfl,r.fo)
    end
    orders=log2.(errors[1:2]./errors[2:3])
    r=last(runs)
    ensure_finite(orders,model,r.steps,r.t_final,r.max_cfl,r.fo)
    return Dict("initial"=>"mode","nx"=>nx,"dx"=>[r.dx for r in runs],"errors"=>errors,"orders"=>orders)
end
function summary_section(r)
    section=Dict{String,Any}("nx"=>length(r.x),"model"=>string(r.model),"advection"=>r.advection,
        "initial"=>string(r.initial),"speed_role"=>r.model==:linear ? "constant advection speed" : "unused; velocity is u, analytic mean is 1")
    for key in (:speed,:diffusivity,:dx,:dt,:steps,:t_final,:requested_safety,:initial_integral,
        :final_integral,:integral_change,:minimum,:maximum,:max_cfl,:fo,:max_stability_number)
        value=getproperty(r,key)
        ensure_finite((value,),r.model,r.steps,r.t_final,r.max_cfl,r.fo)
        section[string(key)]=value
    end
    return section
end
periodic_display(r,u)=(vcat(r.x,2.),vcat(u,u[1]))
const CONDITION_COLORS=("#0072B2","#D55E00","#009E73")
function make_plots(directory,a,b,c,convergence)
    p=plot(periodic_display(c,c.u0)...;label="Initial (t = 0)",color=:gray,linestyle=:dash,linewidth=2,
        xlabel="x (dimensionless)",ylabel=c.model==:linear ? "Temperature u (dimensionless)" : "Velocity u (dimensionless)",
        title="$(c.model), t = $(c.t_final)",size=(800,500),ylims=(0.95,2.05),legend=:outerright)
    for ((r,label),color) in zip(((a,"Advection only"),(b,"Diffusion only"),(c,"Combined")),CONDITION_COLORS)
        plot!(p,periodic_display(r,r.u)...;label,color,linestyle=:solid,linewidth=2)
    end
    savefig(p,joinpath(directory,"comparison.png"))
    conservation_title=c.model==:linear ? "Temperature integral conservation error" : "Velocity integral conservation error"
    p=plot(;xlabel="t (dimensionless)",ylabel="(I - I0) / 1e-14",size=(800,500),
        title="$conservation_title\n$(c.model): I0 = $(round(c.initial_integral;digits=5)) (dimensionless)",
        titlefontsize=12,legend=:outerright)
    for ((r,label),color) in zip(((a,"Advection only"),(b,"Diffusion only"),(c,"Combined")),CONDITION_COLORS)
        plot!(p,r.times,(r.integral_history.-r.initial_integral)./1e-14;label,color,linestyle=:solid,linewidth=2)
    end
    savefig(p,joinpath(directory,"conservation.png"))
    dx=convergence["dx"]; errors=convergence["errors"]
    p=plot(dx,errors;label="Combined ($(c.model))",color=CONDITION_COLORS[3],marker=:circle,linewidth=2,
        xscale=:log10,yscale=:log10,xlabel="dx (dimensionless)",ylabel="Maximum absolute error",
        title="Smooth analytic mode, t = 1",size=(800,500),legend=:topleft)
    plot!(p,dx,errors[1].*dx./dx[1];label="First order",color=:black,linestyle=:dot)
    savefig(p,joinpath(directory,"convergence.png"))
end
const OUTPUT_NAMES=("comparison.png","conservation.png","convergence.png","summary.toml")
function check_output_sizes(staged,output_dir,model)
    target=abspath(joinpath(output_dir,string(model)))
    replaced=Set(joinpath(target,n) for n in OUTPUT_NAMES)
    sizes=[filesize(joinpath(staged,n)) for n in OUTPUT_NAMES]
    all(s->0<s<=5*1024^2,sizes) || error("N04の1ファイル上限5 MiBを超えています")
    retained=0
    if isdir(output_dir)
        for (root,_,files) in walkdir(output_dir), file in files
            path=abspath(joinpath(root,file))
            path in replaced && continue
            size=filesize(path)
            size<=5*1024^2 || error("1ファイル上限5 MiBを超えています: $path")
            retained+=size
        end
    end
    sum(sizes)+retained<=10*1024^2 || error("N04全体の上限10 MiBを超えています")
    # Standard course location: account for other assignments' existing results.
    exercises=dirname(@__DIR__)
    total=sum(sizes)+retained
    if basename(exercises)=="exercises" && isdir(exercises)
        for task in readdir(exercises;join=true)
            results=joinpath(task,"results")
            task==(@__DIR__) && continue
            isdir(results) || continue
            for (root,_,files) in walkdir(results), file in files
                total+=filesize(joinpath(root,file))
            end
        end
    end
    total<=100*1024^2 || error("全課題の出力上限100 MiBを超えています")
end
install_output(source,destination)=cp(source,destination;force=true)
function save_staged(draw,output_dir,model,summary)
    temporary=mktempdir(;cleanup=false)
    preserve_backup=false
    try
        staged=joinpath(temporary,"new"); backup=joinpath(temporary,"backup")
        mkpath(staged); mkpath(backup)
        open(joinpath(staged,"summary.toml"),"w") do io
            TOML.print(io,summary;sorted=true)
        end
        draw(staged)
        check_output_sizes(staged,output_dir,model)
        for n in OUTPUT_NAMES[1:3]
            open(joinpath(staged,n)) do io
                read(io,8)==UInt8[0x89,0x50,0x4e,0x47,0x0d,0x0a,0x1a,0x0a] || error("PNG出力が不正です: $n")
            end
        end
        target=joinpath(output_dir,string(model))
        existing=Set{String}()
        for n in OUTPUT_NAMES
            path=joinpath(target,n)
            ispath(path) && (!isfile(path) || islink(path)) && error("公式出力先は通常ファイルである必要があります: $path")
            if isfile(path)
                cp(path,joinpath(backup,n))
                push!(existing,n)
            end
        end
        mkpath(target)
        touched=String[]
        try
            for n in OUTPUT_NAMES
                push!(touched,n)
                install_output(joinpath(staged,n),joinpath(target,n))
            end
        catch original
            try
                for n in touched
                    destination=joinpath(target,n)
                    if n in existing
                        cp(joinpath(backup,n),destination;force=true)
                    elseif isfile(destination)
                        rm(destination)
                    end
                end
            catch recovery
                preserve_backup=true
                error("出力反映と復元に失敗しました。バックアップ: $backup; 反映: $(sprint(showerror,original)); 復元: $(sprint(showerror,recovery))")
            end
            rethrow(original)
        end
    finally
        preserve_backup || rm(temporary;recursive=true)
    end
end
function write_outputs(output_dir,a,b,c,convergence)
    units=c.model==:linear ? "dimensionless temperature; integral is transported temperature content" :
        "dimensionless Burgers velocity; integral is velocity integral, not heat"
    summary=Dict("course_id"=>"N04","model"=>string(c.model),"boundary"=>"periodic","domain"=>[0.,2.],
        "units"=>units,"initial"=>"pulse","t_final"=>c.t_final,"requested_safety"=>c.requested_safety,
        "advection_only"=>summary_section(a),"diffusion_only"=>summary_section(b),
        "combined"=>summary_section(c),"convergence"=>convergence)
    save_staged(output_dir,c.model,summary) do directory
        make_plots(directory,a,b,c,convergence)
    end
end
