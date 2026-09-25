# 提供: N07セル中心格子と辺順。計算・保存・解析・作図で共有する。
using HDF5, SHA, TOML
const DEFAULT_OUTPUT_DIR=joinpath(@__DIR__,"results")
const OUTPUT_NAMES=("temperature.h5","burgers.h5","summary.toml","temperature_comparison.png","boundary_heat.png","burgers_fields.png","convergence.png","plots.toml")
const OFFICIAL_GRIDS=((40,30),(80,60),(160,120))
const SAVE_TIMES=[0.,.25,.5,.75,1.]
const SIDES=("west","east","south","north")
case_id(nx,ny)="n$(lpad(nx,3,'0'))x$(lpad(ny,3,'0'))"
cell_centers(n,L)=collect(((1:n).-0.5).*(L/n))

require(ok,msg)=ok ? nothing : throw(ArgumentError(msg))
function schedule(time,t_final=1.)
    require(time isa AbstractVector && length(time)>=2 && all(t->t isa Real && isfinite(t),time) && first(time)==0 && last(time)==t_final && all(diff(time).>0),"timeは0始まり・有限・狭義増加・末尾t_finalです")
    Float64.(time)
end
function validate_budget(time,heat,adv,diff)
    require(time isa AbstractVector && length(time)>=2,"timeは2点以上です")
    schedule(time,last(time)); nt=length(time)
    require(heat isa AbstractVector && length(heat)==nt && all(v->v isa Real && isfinite(v),heat),"heatは時刻に対応する有限列です")
    for a in (adv,diff)
        require(a isa AbstractMatrix && size(a)==(4,nt-1) && all(v->v isa Real && isfinite(v),a),"熱輸送積分は有限な(4,nt-1)配列です")
    end
end
const PERIODIC_CASES=("periodic_advection","periodic_diffusion","periodic_combined")
const THERMAL_CASES=(PERIODIC_CASES...,"closed_fixed","closed_insulated","channel_fixed","channel_insulated",("comparison_"*split(c,"_")[2] for c in PERIODIC_CASES)...)
function thermal_config(id)
    b(k,v=0.)=(;kind=k,value=Float64(v))
    allbc(k)=(;west=b(k),east=b(k),south=b(k),north=b(k))
    require(id in THERMAL_CASES,"未知の温度ケース: $id")
    if startswith(id,"periodic_") || startswith(id,"comparison_")
        mode=split(id,"_")[2]; cx,cy=mode=="diffusion" ? (0.,0.) : (1.,.5)
        return (;cx,cy,kappa=mode=="advection" ? 0. : .05,bc=allbc(:periodic),initial="fourier_v1")
    elseif startswith(id,"closed_")
        return (;cx=0.,cy=0.,kappa=.05,bc=allbc(id=="closed_fixed" ? :dirichlet : :insulated),initial="sine_wall_v1")
    end
    wall=id=="channel_fixed" ? :dirichlet : :insulated
    (;cx=1.,cy=0.,kappa=.05,bc=(;west=b(:inflow,1.),east=b(:outflow),south=b(wall),north=b(wall)),initial="uniform_0.2")
end
function exact_temperature(x,y,t,id)
    c=thermal_config(id)
    if c.initial=="fourier_v1"
        X=x.-c.cx*t; Y=y.-c.cy*t; k=c.kappa
        return [1+.2exp(-k*pi^2*t)*sin(pi*xi)+.3exp(-4k*pi^2*t)*cos(2pi*yj)+.1exp(-5k*pi^2*t)*sin(pi*xi+2pi*yj) for xi in X,yj in Y]
    elseif id=="closed_fixed" || (id=="closed_insulated" && t==0)
        return [sin(pi*xi/2)*sin(pi*yj)*exp(-c.kappa*((pi/2)^2+pi^2)*t) for xi in x,yj in y]
    elseif startswith(id,"channel_") && t==0
        return fill(.2,length(x),length(y))
    end
    nothing
end

"""Cole–Hopfの明示微分。X=x-.6t, Y=y+.3t。数値差分で解析値を作らない。"""
function exact_burgers(x,y,t;nu=.05)
    X=x.-.6t; Y=y.+.3t; A=.2exp(-nu*pi^2*t); B=.2exp(-4nu*pi^2*t)
    phi=[1+A*cos(pi*xi)+B*cos(2pi*yj) for xi in X,yj in Y]
    u=[.6+2nu*A*pi*sin(pi*xi)/phi[i,j] for (i,xi) in enumerate(X),(j,yj) in enumerate(Y)]
    v=[-.3+4nu*B*pi*sin(2pi*yj)/phi[i,j] for (i,xi) in enumerate(X),(j,yj) in enumerate(Y)]
    u,v
end

file_sha(path)=bytes2hex(sha256(read(path)))
write_toml(path,doc)=open(io->TOML.print(io,doc;sorted=true),path,"w")
function check_course_capacity(output_dir,new_total)
    # Find the actual student's exercises root, including custom nested result files.
    task=dirname(abspath(output_dir)); exercises=dirname(task)
    basename(exercises)=="exercises" || return
    total=new_total
    for entry in readdir(exercises;join=true)
        entry==task && continue
        results=joinpath(entry,"results")
        isdir(results) || continue
        for (dir,_,files) in walkdir(results),name in files;total+=filesize(joinpath(dir,name));end
    end
    legacy=joinpath(dirname(exercises),"results")
    if isdir(legacy)
        for (dir,_,files) in walkdir(legacy),name in files;total+=filesize(joinpath(dir,name));end
    end
    require(total<=100*1024^2,"全課題の100 MiB上限を超えています")
end
function publish(stage,output_dir,names;copy_file=(a,b)->cp(a,b;force=true),restore_file=(a,b)->cp(a,b;force=true),file_limit=5*1024^2,task_limit=10*1024^2)
    require(all(n->n in OUTPUT_NAMES,names) && length(unique(names))==length(names),"反映できる名前はN07公式出力だけです")
    sizes=[filesize(joinpath(stage,n)) for n in names]
    require(all(sizes.<=file_limit),"出力の1ファイル上限を超えています")
    # Include retained user files in the one N07 total.
    totals=Dict{String,Int}()
    group(n)="N07"
    for (n,size) in zip(names,sizes); totals[group(n)]=get(totals,group(n),0)+size; end
    if isdir(output_dir)
        for (dir,_,files) in walkdir(output_dir), file in files
            path=joinpath(dir,file); rel=relpath(path,output_dir)
            rel in names && continue
            require(filesize(path)<=file_limit,"追加ファイルの1ファイル上限を超えています")
            key=group(rel); totals[key]=get(totals,key,0)+filesize(path)
        end
    end
    require(all(values(totals).<=task_limit),"出力の内容ID合計上限を超えています")
    check_course_capacity(output_dir,sum(values(totals)))
    backup=mktempdir(;cleanup=false); existed=Dict{String,Bool}(); applied=String[]
    try
        for name in names
            target=joinpath(output_dir,name)
            require(!ispath(target) || isfile(target),"出力先が通常ファイルではありません: $target")
            existed[name]=isfile(target)
            if existed[name]
                saved=joinpath(backup,name); mkpath(dirname(saved)); cp(target,saved)
            end
        end
        for name in names
            target=joinpath(output_dir,name); mkpath(dirname(target))
            push!(applied,name) # Include partial writes from a failing copy.
            copy_file(joinpath(stage,name),target)
        end
    catch original
        failures=String[]
        for name in reverse(applied)
            target=joinpath(output_dir,name)
            try
                existed[name] ? restore_file(joinpath(backup,name),target) : rm(target;force=true)
            catch e
                push!(failures,"$target: $e")
            end
        end
        if !isempty(failures)
            error("反映と復元に失敗しました。バックアップ: $backup\n"*join(failures,"\n")*"\n元のエラー: $original")
        end
        rm(backup;recursive=true)
        rethrow()
    end
    rm(backup;recursive=true)
    nothing
end
function staged(action,output_dir,names;kwargs...)
    mktempdir() do temporary
        action(temporary)
        publish(temporary,output_dir,names;kwargs...)
    end
end

const CASE_ATTRIBUTES=("nx","ny","dx","dy","grid_location","cx","cy","kappa","nu","initial_condition","safety","timestep_policy",("$(s)_$(k)" for s in SIDES for k in ("kind","value"))...)
const ROOT_KEYS=("schema_version","task_id","run_id","equation_family","code_sha256","git_head")
const COMMON_DATASETS=("x","y","time","step_time","step_dt","save_step_index")
field_keys(family)=family=="temperature" ? ("temperature",) : ("u","v")
dataset_keys(family)=(COMMON_DATASETS...,field_keys(family)...,(family=="temperature" ? ("advective_heat_integral","diffusive_heat_integral") : ())...)
function source_metadata(run_id,family)
    root=normpath(joinpath(@__DIR__,"..",".."))
    paths=vcat([joinpath("src",n) for n in readdir(joinpath(root,"src")) if endswith(n,".jl")],[relpath(joinpath(@__DIR__,n),root) for n in ("simulate.jl","provided_support.jl")])
    hashes=Dict(n=>file_sha(joinpath(root,n)) for n in paths)
    code=sprint(io->TOML.print(io,hashes;sorted=true))
    head=try readchomp(pipeline(`git -C $root rev-parse HEAD`;stderr=devnull)) catch; "unversioned" end
    Dict{String,Any}("schema_version"=>1,"task_id"=>"N07","run_id"=>run_id,"equation_family"=>family,"code_sha256"=>code,"git_head"=>head)
end
function grids_for(id)
    id in PERIODIC_CASES || id=="periodic_burgers" ? OFFICIAL_GRIDS : ((80,60),)
end
function expected_cases(family)
    Dict(id=>Set(case_id(n...) for n in grids_for(id)) for id in (family=="temperature" ? THERMAL_CASES : ("periodic_burgers",)))
end
function case_metadata(id,nx,ny)
    b=id=="periodic_burgers"
    c=b ? (;cx=.6,cy=-.3,kappa=0.,bc=thermal_config("periodic_advection").bc,initial="cole_hopf_v1") : thermal_config(id)
    doc=Dict{String,Any}("nx"=>nx,"ny"=>ny,"dx"=>2/nx,"dy"=>1/ny,"grid_location"=>"cell_center","cx"=>c.cx,"cy"=>c.cy,"kappa"=>c.kappa,"nu"=>b ? .05 : 0.,"initial_condition"=>c.initial,"safety"=>.8,"timestep_policy"=>b ? "dynamic_old_field" : startswith(id,"comparison_") ? "common_combined" : "case_stable")
    for side in SIDES
        v=getproperty(c.bc,Symbol(side));doc[side*"_kind"]=string(v.kind);doc[side*"_value"]=v.value
    end
    doc
end
read_attribute(o,k)=haskey(attributes(o),k) ? read(attributes(o)[k]) : throw(ArgumentError("必要属性がありません: $k"))
function validate_case(c,id,grid,family)
    nx,ny=c["nx"],c["ny"]
    require(nx isa Integer && ny isa Integer && nx>=3 && ny>=3 && grid==case_id(nx,ny),"$id/$grid: 格子数が不正です")
    expected=case_metadata(id,nx,ny)
    require(all(c[k]==expected[k] for k in CASE_ATTRIBUTES),"$id: 境界・係数・格子配置・刻み方針が不一致です")
    for (key,n,L) in (("x",nx,2.),("y",ny,1.))
        require(c[key] isa Vector{Float64} && length(c[key])==n && all(isfinite,c[key]) && isapprox(c[key],cell_centers(n,L);rtol=1e-13,atol=1e-14),"$id/$key: セル中心座標が不正です")
    end
    time=schedule(c["time"]);require(time==SAVE_TIMES,"公式保存時刻が不一致です"); nt=length(time)
    for key in field_keys(family)
        require(c[key] isa Array{Float64,3} && size(c[key])==(nx,ny,nt) && all(isfinite,c[key]),"$id/$key: shape・有限値が不正です")
    end
    ts,dt,idx=c["step_time"],c["step_dt"],c["save_step_index"]
    require(ts isa Vector{Float64} && dt isa Vector{Float64} && length(ts)==length(dt)>0 && all(isfinite,ts) && all(v->isfinite(v)&&v>0,dt),"$id: step配列が不正です")
    require(first(ts)==0 && all(diff(ts).>0) && all(isapprox.(ts[2:end],ts[1:end-1]+dt[1:end-1];atol=1e-13,rtol=1e-13)),"$id: step時刻とdtが不一致です")
    require(idx isa AbstractVector{<:Integer} && length(idx)==nt && first(idx)==0 && last(idx)==length(dt) && all(diff(idx).>0),"$id: 保存step対応が不正です")
    require(all(isapprox(time[k],ts[idx[k]]+dt[idx[k]];atol=1e-13,rtol=1e-13) for k in 2:nt),"$id: 保存時刻とstepが不一致です")
    if family=="temperature"
        for key in ("advective_heat_integral","diffusive_heat_integral")
            require(c[key] isa Matrix{Float64} && size(c[key])==(4,nt-1) && all(isfinite,c[key]),"$id/$key: 辺・区間・有限値が不正です")
        end
        # Independent metadata stability check; no numerical package is imported for reading.
        ax=any(c[s*"_kind"] in ("dirichlet","inflow") for s in ("west","east")) ? 3 : 2
        ay=any(c[s*"_kind"] in ("dirichlet","inflow") for s in ("south","north")) ? 3 : 2
        rate=abs(c["cx"])/c["dx"]+abs(c["cy"])/c["dy"]+c["kappa"]*(ax/c["dx"]^2+ay/c["dy"]^2)
        startswith(id,"comparison_") && (rate=1/c["dx"]+.5/c["dy"]+.1*(1/c["dx"]^2+1/c["dy"]^2))
        require(all(dt.*rate .<=.8+32eps()),"$id: 安定刻み超過です")
    end
end
function write_fields(path,cases,meta;official=true)
    family=meta["equation_family"]
    h5open(path,"w") do h
        for (k,v) in meta;attributes(h)[k]=v;end
        cg=create_group(h,"cases")
        for (id,grids) in sort(collect(cases);by=first)
            ig=create_group(cg,id)
            for (grid,c) in sort(collect(grids);by=first)
                g=create_group(ig,grid)
                for k in CASE_ATTRIBUTES;attributes(g)[k]=c[k];end
                for k in dataset_keys(family)
                    g[k]=c[k]
                    k!="save_step_index" && (attributes(g[k])["units"]="1")
                    k in field_keys(family) && (attributes(g[k])["axis_order"]="time,y,x")
                    if endswith(k,"heat_integral")
                        attributes(g[k])["axis_order"]="interval,side";attributes(g[k])["side_order"]="west,east,south,north";attributes(g[k])["positive_direction"]="inward"
                    end
                end
            end
        end
    end
    read_fields(path;official)
end
function read_fields(path;official=true)
    isfile(path) || error("HDF5欠落: $(path)。simulate.jlを実行してください")
    try
        h5open(path,"r") do h
            meta=Dict(k=>read_attribute(h,k) for k in ROOT_KEYS); family=meta["equation_family"]
            require(meta["schema_version"]==1 && meta["task_id"]=="N07" && family in ("temperature","burgers"),"root schema不一致")
            require(meta["run_id"] isa String && occursin(r"^[0-9a-f]{32}$",meta["run_id"]),"run_id不正")
            hashes=TOML.parse(meta["code_sha256"])
            require(!isempty(hashes) && all(v->v isa String && occursin(r"^[0-9a-f]{64}$",v),values(hashes)) && meta["git_head"] isa String && !isempty(meta["git_head"]),"コード来歴が不正です")
            require(haskey(h,"cases"),"cases欠落");cases=Dict{String,Any}()
            for id in keys(h["cases"])
                require(id in keys(expected_cases(family)),"未知ケースです: $id")
                grids=Dict{String,Any}()
                for grid in keys(h["cases/$id"])
                    g=h["cases/$id/$grid"]; c=Dict{String,Any}(k=>read_attribute(g,k) for k in CASE_ATTRIBUTES)
                    for k in dataset_keys(family)
                        require(haskey(g,k),"$id/$grid/$(k)欠落");c[k]=read(g[k])
                        k!="save_step_index" && require(read_attribute(g[k],"units")=="1","units不一致")
                        k in field_keys(family) && require(read_attribute(g[k],"axis_order")=="time,y,x","場の軸不一致")
                        if endswith(k,"heat_integral")
                            require(read_attribute(g[k],"axis_order")=="interval,side" && read_attribute(g[k],"side_order")=="west,east,south,north" && read_attribute(g[k],"positive_direction")=="inward","熱輸送の軸・辺順・符号が不一致")
                        end
                    end
                    validate_case(c,id,grid,family);grids[grid]=c
                end
                cases[id]=grids
            end
            if official
                expected=expected_cases(family)
                require(Set(keys(cases))==Set(keys(expected)) && all(Set(keys(cases[id]))==expected[id] for id in keys(expected)),"公式ケース・格子が不足しています")
            end
            (;metadata=meta,cases)
        end
    catch e
        error("HDF5読取り失敗 $path: $(sprint(showerror,e))")
    end
end
function read_pair(dir)
    t=read_fields(joinpath(dir,"temperature.h5")); b=read_fields(joinpath(dir,"burgers.h5"))
    require(t.metadata["equation_family"]=="temperature" && b.metadata["equation_family"]=="burgers","HDF5の方程式系が不一致")
    for k in ("run_id","code_sha256","git_head");require(t.metadata[k]==b.metadata[k],"HDF5の計算由来が混在しています: $k");end
    (;temperature=t,burgers=b)
end
function input_hashes(dir)
    Dict(n=>file_sha(joinpath(dir,n)) for n in ("temperature.h5","burgers.h5"))
end

function diagnostics(pair,budget)
    cases=Dict{String,Any}();convergence=Dict{String,Any}()
    for (family,data) in (("temperature",pair.temperature),("burgers",pair.burgers))
        for (id,grids) in data.cases
            ds=Dict{String,Any}()
            for (grid,c) in grids
                d=Dict{String,Any}("conditions"=>Dict(k=>c[k] for k in CASE_ATTRIBUTES),"time"=>c["time"],"step_count"=>length(c["step_dt"]),"dt_min"=>minimum(c["step_dt"]),"dt_max"=>maximum(c["step_dt"]))
                if family=="temperature"
                    U=c["temperature"];d["heat"]=[c["dx"]*c["dy"]*sum(@view U[:,:,k]) for k in eachindex(c["time"])]
                    for (key,fun) in (("minimum",minimum),("maximum",maximum));d[key]=[fun(@view U[:,:,k]) for k in eachindex(c["time"])];end
                    adv=c["advective_heat_integral"];diff=c["diffusive_heat_integral"]
                    b=budget(c["time"],d["heat"],adv,diff)
                    for (k,v) in pairs(b);require(all(isfinite,v),"$id: 非有限な熱収支診断です");d[string(k)]=v;end
                    d["advective_heat_integral"]=Dict(side=>collect(adv[j,:]) for (j,side) in enumerate(SIDES))
                    d["diffusive_heat_integral"]=Dict(side=>collect(diff[j,:]) for (j,side) in enumerate(SIDES))
                    if id in PERIODIC_CASES || startswith(id,"comparison_") || id=="closed_fixed"
                        d["l2_error"]=[sqrt(sum(abs2,U[:,:,k]-exact_temperature(c["x"],c["y"],t,id))/(c["nx"]*c["ny"])) for (k,t) in enumerate(c["time"])]
                    end
                else
                    for component in ("u","v")
                        U=c[component];d[component*"_minimum"]=[minimum(@view U[:,:,k]) for k in eachindex(c["time"])];d[component*"_maximum"]=[maximum(@view U[:,:,k]) for k in eachindex(c["time"])];d[component*"_l2_error"]=Float64[]
                    end
                    for (k,t) in enumerate(c["time"])
                        ue,ve=exact_burgers(c["x"],c["y"],t)
                        for (key,exact) in (("u",ue),("v",ve));push!(d[key*"_l2_error"],sqrt(sum(abs2,c[key][:,:,k]-exact)/length(exact)));end
                    end
                    d["l2_error"]=hypot.(d["u_l2_error"],d["v_l2_error"])
                end
                ds[grid]=d
            end
            cases[id]=ds
            if id in PERIODIC_CASES || id=="periodic_burgers"
                ids=sort(collect(keys(ds));by=g->ds[g]["conditions"]["nx"])
                errors=[last(ds[g]["l2_error"]) for g in ids]
                require(all(v->isfinite(v)&&v>0,errors),"最終誤差が不正です")
                convergence[id]=Dict("grid_ids"=>ids,"errors"=>errors,"orders"=>log2.(errors[1:end-1]./errors[2:end]))
            end
        end
    end
    Dict{String,Any}("schema_version"=>1,"task_id"=>"N07","run_id"=>pair.temperature.metadata["run_id"],"diagnostics_complete"=>true,"provenance"=>pair.temperature.metadata,"cases"=>cases,"convergence"=>convergence)
end
function read_summary(path,input_dir)
    require(isfile(path),"summary欠落。analyze.jlを実行してください")
    s=TOML.parsefile(path);data=read_pair(input_dir)
    require(get(s,"schema_version",0)==1 && get(s,"task_id","")=="N07" && get(s,"source_sha256",nothing)==input_hashes(input_dir) && get(s,"run_id",nothing)==data.temperature.metadata["run_id"],"summaryの出自がHDF5と不一致です。analyze.jlを再実行してください")
    require(get(s,"diagnostics_complete",false)===true,"熱収支解析が未完了です")
    require(haskey(s,"cases") && Set(keys(s["cases"]))==union(Set(keys(data.temperature.cases)),Set(keys(data.burgers.cases))),"summaryのケースが不足しています")
    s
end
function check_complete(output_dir=DEFAULT_OUTPUT_DIR)
    pair=read_pair(output_dir);s=read_summary(joinpath(output_dir,"summary.toml"),output_dir)
    for (id,conv) in s["convergence"]
        expected=grids_for(id);require(length(conv["errors"])==length(expected)==3,"3格子が必要です")
        errors=conv["errors"];p=conv["orders"];lo,hi=id=="periodic_diffusion" ? (1.8,2.2) : (.8,1.2)
        require(all(diff(errors).<0) && length(p)==2 && all(lo .<= p .<= hi) && isapprox(p,log2.(errors[1:2]./errors[2:3]);rtol=1e-12),"$id: 収束基準を満たしません")
    end
    require(Set(keys(s["convergence"]))==Set((PERIODIC_CASES...,"periodic_burgers")),"収束条件が不足しています")
    # Verify student diagnostics against persisted inputs, without rerunning a solver.
    for (id,grids) in pair.temperature.cases, (grid,c) in grids
        d=s["cases"][id][grid];q=[c["dx"]*c["dy"]*sum(@view c["temperature"][:,:,k]) for k in eachindex(c["time"])]
        transport=sum(abs,c["advective_heat_integral"])+sum(abs,c["diffusive_heat_integral"])
        tol=1e-12*max(1,abs(first(q)),transport)
        require(isapprox(d["heat"],q;rtol=1e-13,atol=1e-14),"$id: 熱量が保存場と不一致です")
        require(length(d["residual"])==length(q) && maximum(abs,d["residual"])<=tol && isapprox(q.-first(q),d["cumulative_input"]+d["residual"];atol=tol,rtol=0),"$id: 熱収支が不一致です")
        require(first(d["cumulative_input"])==0 && isapprox(diff(d["cumulative_input"]),d["net_input"];rtol=1e-12,atol=tol),"$id: 区間と累積収支が不一致です")
        for k in eachindex(d["net_input"])
            require(isapprox(d["net_input"][k],sum(c["advective_heat_integral"][:,k])+sum(c["diffusive_heat_integral"][:,k]);atol=tol,rtol=0),"$id: 境界輸送と不一致です")
        end
    end
    # Recompute the analytic error from stored fields; a small heat residual alone is insufficient.
    for id in (PERIODIC_CASES...,"periodic_burgers")
        data=id=="periodic_burgers" ? pair.burgers : pair.temperature
        conv=s["convergence"][id]
        for (k,grid) in enumerate(conv["grid_ids"])
            c=data.cases[id][grid]
            if id=="periodic_burgers"
                u,v=exact_burgers(c["x"],c["y"],1.);e=sqrt((sum(abs2,c["u"][:,:,end]-u)+sum(abs2,c["v"][:,:,end]-v))/length(u))
            else
                u=exact_temperature(c["x"],c["y"],1.,id);e=sqrt(sum(abs2,c["temperature"][:,:,end]-u)/length(u))
            end
            require(isapprox(conv["errors"][k],e;rtol=1e-12),"$id: 解析解誤差が保存場と不一致です")
        end
    end
    path=joinpath(output_dir,"plots.toml");require(isfile(path),"plots.toml欠落")
    p=TOML.parsefile(path)
    require(p["source_sha256"]==input_hashes(output_dir) && p["source_summary_sha256"]==file_sha(joinpath(output_dir,"summary.toml")),"図の出自が不一致です。plot.jlを再実行してください")
    require(all(n->isfile(joinpath(output_dir,n)),OUTPUT_NAMES),"公式8出力が不足しています")
    require(all(file_sha(joinpath(output_dir,n))==p["figure_sha256"][n] for n in keys(p["figure_sha256"])) && Set(keys(p["figure_sha256"]))==Set(n for n in OUTPUT_NAMES if endswith(n,".png")),"図が記録と不一致です")
    true
end
