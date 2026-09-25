# 提供: 保存・検証・解析解・出力反映。数値更新はsrc/、追加診断はanalyze.jlで実装する。
using HDF5, SHA, TOML
const DEFAULT_OUTPUT_DIR=joinpath(@__DIR__,"results","N06")
const FIELD_NAME="fields.h5"
const OUTPUT_NAMES=("fields.h5","summary.toml","fields.png","diagnostics.png","convergence.png","plots.toml")
const OFFICIAL_GRIDS=((40,30),(80,60),(160,120))
const SAVE_TIMES=[0.,.25,.5,.75,1.]
const ROOT_ATTRIBUTES=Dict{String,Any}("schema_version"=>1,"task_id"=>"N06","equation"=>"linear_advection_2d",
    "initial_condition"=>"smooth_periodic_v1","boundary_x"=>"periodic","boundary_y"=>"periodic",
    "cx"=>1.,"cy"=>.5,"Lx"=>2.,"Ly"=>1.,"safety"=>.8,"t_final"=>1.)
file_sha(path)=bytes2hex(sha256(read(path)))
case_id(nx,ny)="n$(lpad(nx,3,'0'))x$(lpad(ny,3,'0'))"
u0(x,y)=1+.2sin(pi*x)+.3cos(2pi*y)+.1sin(pi*x+2pi*y)
exact_field(x,y,t;cx=1.,cy=.5)=[u0(mod(xi-cx*t,2.),mod(yj-cy*t,1.)) for xi in x,yj in y]
require(ok,message)=ok ? nothing : throw(ArgumentError(message))
function schedule(times,t_final)
    require(times isa AbstractVector && length(times)>=2,"time: 2時刻以上が必要です")
    require(all(t->t isa Real && isfinite(t),times) && first(times)==0 && last(times)==t_final && all(diff(times).>0),"time: 有限・狭義増加・先頭0・末尾t_finalが必要です")
    Float64.(times)
end
function write_toml(path,doc)
    open(io->TOML.print(io,doc;sorted=true),path,"w")
end
"""全出力を生成してから反映する。反映失敗時は既存ファイルを復元する。"""
function publish(stage,output_dir,names;copy_file=(a,b)->cp(a,b;force=true),restore_file=(a,b)->cp(a,b;force=true),file_limit=5*1024^2,task_limit=10*1024^2)
    sizes=[filesize(joinpath(stage,n)) for n in names]
    require(all(sizes.<=file_limit),"出力の1ファイル上限を超えています")
    # Count retained files as well; combined run uses N05/ and N06/ as separate content IDs.
    totals=Dict{String,Int}()
    group(n)=first(splitpath(n)) in ("N05","N06") ? first(splitpath(n)) : "N06"
    for (n,size) in zip(names,sizes); totals[group(n)]=get(totals,group(n),0)+size; end
    if isdir(output_dir)
        for (dir,_,files) in walkdir(output_dir), file in files
            path=joinpath(dir,file); rel=relpath(path,output_dir)
            rel in names && continue
            key=group(rel); totals[key]=get(totals,key,0)+filesize(path)
        end
    end
    require(all(values(totals).<=task_limit),"出力の内容ID合計上限を超えています")
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
function write_fields(path,cases;metadata=ROOT_ATTRIBUTES)
    h5open(path,"w") do h
        for (k,v) in metadata; attributes(h)[k]=v; end
        group=create_group(h,"cases")
        for (id,c) in sort(collect(cases);by=first)
            g=create_group(group,id)
            for key in ("nx","ny","dx","dy","dt_max"); attributes(g)[key]=c[key]; end
            for key in ("x","y","time","u","segment_dt","segment_steps")
                g[key]=c[key]
                key!="segment_steps" && (attributes(g[key])["units"]="1")
            end
            attributes(g["u"])["axis_order"]="time,y,x"
        end
    end
    # Validate real persisted values before publishing.
    read_fields(path)
end
read_attribute(obj,key)=haskey(attributes(obj),key) ? read(attributes(obj)[key]) : error("必要属性がありません: $key")
function validate_case(c,id,meta)
    nx,ny=c["nx"],c["ny"]
    require(nx isa Integer && ny isa Integer && nx>=3 && ny>=3,"$id: nx,nyは3以上の整数です")
    require(id==case_id(nx,ny),"$id: ケース名と格子数が不一致です")
    dx,dy=c["dx"],c["dy"]
    require(all(v->v isa Real && isfinite(v) && v>0,(dx,dy,c["dt_max"])),"$id: 格子幅・刻み上限が不正です")
    require(isapprox(dx,meta["Lx"]/nx;rtol=1e-13) && isapprox(dy,meta["Ly"]/ny;rtol=1e-13),"$id: 格子幅が領域と不一致です")
    require(c["time"] isa Vector{Float64},"$id/time: Float64配列が必要です")
    time=schedule(c["time"],meta["t_final"]); nt=length(time)
    for (key,n,h) in (("x",nx,dx),("y",ny,dy))
        v=c[key]
        require(v isa Vector{Float64} && length(v)==n && all(isfinite,v) && all(isapprox.(v,(0:n-1).*h;rtol=1e-13,atol=1e-14)),"$id/$key: 座標・周期終点・形状が不一致です")
    end
    U=c["u"]
    require(U isa Array{Float64,3} && size(U)==(nx,ny,nt) && all(isfinite,U),"$id/u: 軸・形状または有限値が不正です")
    dts,steps=c["segment_dt"],c["segment_steps"]
    require(dts isa Vector{Float64} && steps isa AbstractVector{<:Integer} && length(dts)==length(steps)==nt-1,"$id: 保存区間の配列形状が不正です")
    require(all(v->isfinite(v)&&v>0,dts) && all(>(0),steps),"$id: 実効刻み・ステップ数が不正です")
    require(all(isapprox.(dts.*steps,diff(time);rtol=1e-13,atol=1e-14)),"$id: 刻み×ステップ数と保存時刻差が不一致です")
    rate=meta["cx"]/dx+meta["cy"]/dy
    require(isapprox(c["dt_max"],meta["safety"]/rate;rtol=1e-13),"$id: dt_maxと条件が不一致です")
    require(all(dts.*rate .<=1+32eps()) && all(dts .<= c["dt_max"]*(1+32eps())),"$id: 合成CFLまたは刻み上限を超えています")
    nothing
end
function read_fields(path)
    isfile(path) || error("入力HDF5欠落: $(path)。simulate.jlを実行してください。")
    try
        h5open(path,"r") do h
            meta=Dict(k=>read_attribute(h,k) for k in keys(ROOT_ATTRIBUTES))
            for key in ("schema_version","task_id","equation","initial_condition","boundary_x","boundary_y","Lx","Ly")
                require(meta[key]==ROOT_ATTRIBUTES[key],"root/$key: schema不一致です")
            end
            require(all(v->v isa Real && isfinite(v) && v>=0,(meta["cx"],meta["cy"])) && meta["cx"]+meta["cy"]>0,"root: 速度が不正です")
            require(meta["safety"] isa Real && 0<meta["safety"]<=1 && meta["t_final"] isa Real && isfinite(meta["t_final"]) && meta["t_final"]>0,"root: safety/t_finalが不正です")
            require(haskey(h,"cases") && !isempty(keys(h["cases"])),"casesがありません")
            cases=Dict{String,Any}()
            for id in keys(h["cases"])
                g=h["cases/$id"]; c=Dict{String,Any}(k=>read_attribute(g,k) for k in ("nx","ny","dx","dy","dt_max"))
                for key in ("x","y","time","u","segment_dt","segment_steps")
                    require(haskey(g,key),"$id/$key: dataset欠落")
                    c[key]=read(g[key])
                    key!="segment_steps" && require(read_attribute(g[key],"units")=="1","$id/$key: units不一致")
                end
                require(read_attribute(g["u"],"axis_order")=="time,y,x","$id/u: axis_order不一致")
                validate_case(c,id,meta); cases[id]=c
            end
            (;metadata=meta,cases)
        end
    catch e
        error("HDF5読取り失敗: $(path): $(sprint(showerror,e))")
    end
end
function diagnostics(data,variance)
    cases=Dict{String,Any}(); complete=true
    for (id,c) in data.cases
        d=Dict{String,Any}(k=>c[k] for k in ("nx","ny","dx","dy","dt_max","time","segment_dt","segment_steps"))
        for key in ("mass","mean","minimum","maximum","l2_error","linf_error","variance"); d[key]=Float64[]; end
        for (k,t) in enumerate(c["time"])
            u=c["u"][:,:,k]; exact=exact_field(c["x"],c["y"],t;cx=data.metadata["cx"],cy=data.metadata["cy"])
            for (key,value) in (("mass",c["dx"]*c["dy"]*sum(u)),("mean",sum(u)/length(u)),("minimum",minimum(u)),("maximum",maximum(u)),("l2_error",sqrt(sum(abs2,u-exact)/length(u))),("linf_error",maximum(abs,u-exact)))
                require(isfinite(value),"$id/$key: 非有限な診断値です"); push!(d[key],value)
            end
            v=variance(u)
            if isnothing(v)
                complete=false
            else
                require(v isa Float64 && isfinite(v) && v>=0,"$id/variance: 有限・非負のFloat64または未完成のnothingが必要です")
                push!(d["variance"],v)
            end
        end
        cases[id]=d
    end
    if !complete
        for d in values(cases); delete!(d,"variance"); end
    end
    ids=sort(collect(keys(cases));by=id->cases[id]["nx"])
    errors=[last(cases[id]["l2_error"]) for id in ids]
    require(all(>(0),errors),"最終誤差は正である必要があります（初期時刻で収束を評価しません）")
    orders=log2.(errors[1:end-1]./errors[2:end])
    Dict{String,Any}("schema_version"=>1,"diagnostics_complete"=>complete,"conditions"=>data.metadata,"cases"=>cases,
        "convergence"=>Dict("case_ids"=>ids,"errors"=>errors,"orders"=>orders))
end
function read_summary(path,fields)
    isfile(path) || error("解析結果欠落: $(path)。analyze.jlを実行してください。")
    s=try TOML.parsefile(path) catch e; error("解析結果読取り失敗: $(path): $e") end
    require(get(s,"schema_version",0)==1 && get(s,"source_fields_sha256",nothing)==file_sha(fields),"解析結果が現在のHDF5と不一致です: $(path)。analyze.jlを再実行してください。")
    s
end
function check_complete(output_dir=DEFAULT_OUTPUT_DIR)
    fields=joinpath(output_dir,"fields.h5"); data=read_fields(fields)
    summary=joinpath(output_dir,"summary.toml"); s=read_summary(summary,fields)
    require(s["diagnostics_complete"],"spatial_varianceが未実装です")
    require(data.metadata==ROOT_ATTRIBUTES,"公式計算条件と不一致です")
    require(Set(keys(data.cases))==Set(case_id(n...) for n in OFFICIAL_GRIDS),"公式3格子が必要です")
    # Fixed official reference values validate the student diagnostic without supplying its algorithm.
    reference_variance=Dict(
        "n040x030"=>[.07,.06376010565899816,.05810024176767521,.05296472555517874,.04830334553165218],
        "n080x060"=>[.07,.06680655185569011,.0637652648362681,.06086864410038924,.05810957283257867],
        "n160x120"=>[.07,.06838411939118463,.06680722926541642,.06526835656917583,.06376655312125709],
    )
    expected=diagnostics(data,u->nothing)
    for id in keys(expected["cases"]), key in ("time","mass","mean","minimum","maximum","l2_error","linf_error")
        require(haskey(s["cases"][id],key) && isapprox(s["cases"][id][key],expected["cases"][id][key];rtol=1e-12,atol=1e-13),"$id/$key: 解析結果が保存場と不一致です")
    end
    for (id,reference) in reference_variance
        require(haskey(s["cases"][id],"variance") && isapprox(s["cases"][id]["variance"],reference;rtol=1e-12,atol=1e-13),"$id/variance: 公式条件の分散と不一致です")
    end
    for c in values(s["cases"])
        require(c["time"]==SAVE_TIMES,"公式保存時刻が必要です")
        require(maximum(abs.(c["mass"].-first(c["mass"])))<=1e-12*max(1.,abs(first(c["mass"]))),"保存量が変化しています")
    end
    conv=s["convergence"]
    require(conv["case_ids"]==expected["convergence"]["case_ids"] && isapprox(conv["errors"],expected["convergence"]["errors"];rtol=1e-12) && isapprox(conv["orders"],expected["convergence"]["orders"];rtol=1e-12),"収束診断が保存場と不一致です")
    require(all(diff(conv["errors"]).<0) && all(.8 .<=conv["orders"].<=1.2),"3格子の約1次収束が必要です")
    plots=joinpath(output_dir,"plots.toml")
    isfile(plots) || error("plots.toml欠落。plot.jlを実行してください。")
    p=TOML.parsefile(plots)
    require(p["source_fields_sha256"]==file_sha(fields) && p["source_summary_sha256"]==file_sha(summary),"図の出自が古くなっています。plot.jlを再実行してください。")
    require(all(n->isfile(joinpath(output_dir,n)),OUTPUT_NAMES),"公式出力が不足しています")
    true
end
