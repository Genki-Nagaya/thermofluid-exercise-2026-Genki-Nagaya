module N05Regression
using TOML, SHA
const ROOT=normpath(joinpath(@__DIR__,"..",".."))
const DEFAULT_OUTPUT_DIR=joinpath(@__DIR__,"results","N05")
const DIRECTORIES=("N01_linear_advection","N02_nonlinear_advection","N03_diffusion","N04_advection_diffusion")
sha(path)=bytes2hex(sha256(read(path)))
plain(x::NamedTuple)=Dict(string(k)=>plain(v) for (k,v) in pairs(x))
plain(x::AbstractArray)=map(plain,x)
plain(x::Symbol)=string(x)
plain(x)=x
function code_hashes(root)
    paths=String[]
    for dir in vcat([joinpath(root,"src")],[joinpath(root,"exercises",d) for d in DIRECTORIES])
        for (base,_,names) in walkdir(dir), name in names
            endswith(name,".jl") && push!(paths,joinpath(base,name))
        end
    end
    Dict(relpath(p,root)=>sha(p) for p in sort(paths))
end
function capture(root)
    sandbox=Module(gensym(:Baseline))
    for dir in DIRECTORIES
        Base.include(sandbox,joinpath(root,"exercises",dir,"run.jl"))
    end
    Base.invokelatest(capture_loaded,root,sandbox)
end
function capture_loaded(root,sandbox)
    a,b,c,d=(getfield(sandbox,n) for n in (:N01LinearAdvection,:N02NonlinearAdvection,:N03Diffusion,:N04AdvectionDiffusion))
    # Called in a new process; invokelatest bridges included module definitions.
    call(f;kwargs...)=Base.invokelatest(f;kwargs...)
    model=d.SELECTED_MODEL
    results=Dict{String,Any}()
    for scheme in (:upwind,:centered)
        results["N01_$scheme"]=plain(call(a.simulate;scheme))
    end
    for boundary in (:fixed,:periodic)
        results["N02_$boundary"]=plain(call(b.simulate;boundary,nx=boundary==:fixed ? 81 : 80))
    end
    results["N02_high_cfl"]=plain(call(b.simulate;boundary=:periodic,nx=40,cfl=1.1,t_final=.02))
    for boundary in (:fixed,:insulated), initial in (:pulse,:mode)
        results["N03_$(boundary)_$initial"]=plain(call(c.simulate;boundary,initial))
    end
    results["N03_high_fo"]=plain(call(c.simulate;fo=.6,t_final=.01))
    combined=call(d.simulate;model)
    results["N04_combined"]=plain(combined)
    results["N04_advection"]=plain(call(d.simulate;model,diffusivity=0.,dt=combined.dt))
    results["N04_diffusion"]=plain(call(d.simulate;model,advection=false,dt=combined.dt))
    results["N04_mode"]=plain(call(d.simulate;model,initial=:mode))
    head=try strip(read(`git -C $root rev-parse HEAD`,String)) catch; "unavailable" end
    Dict("schema_version"=>1,"selected_model"=>string(model),"git_head"=>head,"code_sha256"=>code_hashes(root),"results"=>results)
end
function fresh_capture(root)
    mktempdir() do tmp
        path=joinpath(tmp,"capture.toml")
        command=`$(Base.julia_cmd()) --startup-file=no --project=$root $(@__FILE__) _capture $root $path`
        run(command)
        TOML.parsefile(path)
    end
end
function compare!(differences,a,b,path="results")
    if a isa AbstractDict
        b isa AbstractDict && Set(keys(a))==Set(keys(b)) || error("N05 キー不一致: $path")
        for k in sort(collect(keys(a))); compare!(differences,a[k],b[k],"$path.$k"); end
    elseif a isa AbstractArray
        b isa AbstractArray && size(a)==size(b) || error("N05 形状不一致: $path")
        for i in eachindex(a); compare!(differences,a[i],b[i],"$path[$i]"); end
    elseif a isa AbstractFloat
        b isa AbstractFloat && isfinite(a) && isfinite(b) && isapprox(a,b;rtol=1e-12,atol=1e-13) || error("N05 数値不一致: $path ($a != $b)")
        differences[path]=abs(a-b)
    else
        typeof(a)==typeof(b) && a==b || error("N05 値不一致: $path ($a != $b)")
    end
end
function write_toml(path,data)
    mkpath(dirname(path))
    temporary,io=mktemp(dirname(path))
    try
        TOML.print(io,data;sorted=true); close(io)
        mv(temporary,path;force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end
function baseline(;root=ROOT,output_dir=DEFAULT_OUTPUT_DIR)
    path=joinpath(output_dir,"baseline.toml")
    ispath(path) && error("N05 baselineは上書きできません: $path")
    data=fresh_capture(root)
    ispath(path) && error("N05 baselineが既に存在します: $path")
    write_toml(path,data)
    path
end
function verify(;root=ROOT,baseline_path=joinpath(DEFAULT_OUTPUT_DIR,"baseline.toml"),output_dir=DEFAULT_OUTPUT_DIR)
    isfile(baseline_path) || error("N05 baseline欠落: $(baseline_path)。抽出前commitの隔離コピーでN05.jl baselineを実行してください。")
    before=try TOML.parsefile(baseline_path) catch e; error("N05 baseline破損: $baseline_path: $e") end
    get(before,"schema_version",0)==1 || error("N05 baseline schema不一致: $baseline_path")
    now=fresh_capture(root)
    get(before,"selected_model",nothing)==now["selected_model"] || error("N05 N04選択変更: $baseline_path")
    differences=Dict{String,Float64}()
    compare!(differences,before["results"],now["results"])
    report=Dict("schema_version"=>1,"baseline_sha256"=>sha(baseline_path),"code_sha256"=>now["code_sha256"],"git_head"=>now["git_head"],"selected_model"=>now["selected_model"],"differences"=>differences)
    write_toml(joinpath(output_dir,"regression.toml"),report)
    report
end
function main(args=ARGS)
    isempty(args) && return verify()
    args==["baseline"] && return baseline()
    args==["verify"] && return verify()
    if length(args)==3 && args[1]=="_capture"
        return write_toml(args[3],capture(args[2]))
    end
    error("使い方: N05.jl [baseline|verify]")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
