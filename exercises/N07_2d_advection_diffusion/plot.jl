module N07Plots
include("provided_support.jl")
using Plots
const COLORMAP=:viridis
const COLOR_LIMITS=(0.0,1.6)
const DISPLAY_GRID="n080x060"
function heatpanel(c,a,title;clims=COLOR_LIMITS,color=COLORMAP)
    heatmap(c["x"],c["y"],permutedims(a);title,xlabel="x",ylabel="y",color,clims,xlims=(0.,2.),ylims=(0.,1.),aspect_ratio=:equal,titlefontsize=10)
end
function draw(stage,pair,s;colormap=COLORMAP,color_limits=COLOR_LIMITS)
    require(length(color_limits)==2 && all(isfinite,color_limits) && color_limits[1]<color_limits[2],"色範囲が不正です")
    temp=pair.temperature.cases;grid=DISPLAY_GRID; panels=[]
    for id in ("comparison_advection","comparison_diffusion","comparison_combined")
        c=temp[id][grid];push!(panels,heatpanel(c,c["temperature"][:,:,end],replace(id,"comparison_"=>"")*" T, t=1";clims=color_limits,color=colormap))
    end
    savefig(plot(panels...;layout=(3,1),size=(850,830),margin=5Plots.mm),joinpath(stage,"temperature_comparison.png"))
    panels=[]
    for id in ("closed_fixed","closed_insulated","channel_fixed","channel_insulated")
        c=temp[id][grid];push!(panels,heatpanel(c,c["temperature"][:,:,end],replace(id,"_"=>" ")*" T";clims=color_limits,color=colormap))
    end
    rates=plot(;title="Boundary transport, inward positive",ylabel="Integrated heat (0 to 1)",legend=:outerright,titlefontsize=10)
    labels=["fixed "*s for s in SIDES];append!(labels,["insulated "*s for s in SIDES])
    av=Float64[];dv=Float64[]
    for id in ("channel_fixed","channel_insulated"),side in SIDES
        d=s["cases"][id][grid];push!(av,sum(d["advective_heat_integral"][side]));push!(dv,sum(d["diffusive_heat_integral"][side]))
    end
    scatter!(rates,1:8,av;label="advection",marker=:circle);scatter!(rates,1:8,dv;label="diffusion",marker=:diamond);plot!(rates;xticks=(1:8,labels),xrotation=35,tickfontsize=7)
    budget=plot(;title="Heat change and accumulated inflow",xlabel="time",ylabel="Heat",legend=:topleft,titlefontsize=10)
    residual=plot(;title="Heat-budget residual",xlabel="time",ylabel="Residual / 1e-15",titlefontsize=10)
    for id in ("channel_fixed","channel_insulated")
        d=s["cases"][id][grid];label=replace(id,"channel_"=>"")
        plot!(budget,d["time"],d["heat"].-first(d["heat"]);label=label*" delta Q",marker=:circle)
        plot!(budget,d["time"],d["cumulative_input"];label=label*" inflow",linestyle=:dash,linewidth=2)
        plot!(residual,d["time"],d["residual"]./1e-15;label,marker=:circle)
    end
    push!(panels,rates,budget,residual)
    # Deliberately separated heat, mechanism and residual panels.
    savefig(plot(panels...;layout=@layout([a b; c d; e{0.4h}; f g]),size=(1250,1300),margin=7Plots.mm),joinpath(stage,"boundary_heat.png"))
    c=pair.burgers.cases["periodic_burgers"][grid];ue,ve=exact_burgers(c["x"],c["y"],1.);panels=[]
    errors=Dict{String,Any}();field_limits=Dict{String,Any}()
    for (key,exact) in (("u",ue),("v",ve))
        numeric=c[key][:,:,end];amp=max(maximum(abs,numeric-exact),eps());bounds=extrema(vcat(vec(numeric),vec(exact)));errors[key]=[-amp,amp];field_limits[key]=collect(bounds)
        push!(panels,heatpanel(c,numeric,"$key numerical, t=1";clims=bounds,color=colormap),heatpanel(c,exact,"$key exact, t=1";clims=bounds,color=colormap),heatpanel(c,numeric-exact,"$key error";clims=(-amp,amp),color=:RdBu))
    end
    savefig(plot(panels...;layout=(2,3),size=(1350,620),margin=5Plots.mm),joinpath(stage,"burgers_fields.png"))
    panels=[]
    for id in (PERIODIC_CASES...,"periodic_burgers")
        conv=s["convergence"][id];h=[s["cases"][id][g]["conditions"]["dx"] for g in conv["grid_ids"]];err=conv["errors"];order=id=="periodic_diffusion" ? 2 : 1
        p=plot(h,err;xscale=:log10,yscale=:log10,title=replace(id,"periodic_"=>""),xlabel="dx (dy refined together)",ylabel="Final RMS L2",label="numerical",marker=:circle,legend=:topleft,titlefontsize=10)
        plot!(p,h,first(err).*(h./first(h)).^order;label="order $order",linestyle=:dash);push!(panels,p)
    end
    savefig(plot(panels...;layout=(2,2),size=(1000,760),margin=6Plots.mm),joinpath(stage,"convergence.png"))
    (;errors,field_limits)
end
function main(;input_dir=DEFAULT_OUTPUT_DIR,summary_path=joinpath(input_dir,"summary.toml"),output_dir=DEFAULT_OUTPUT_DIR,colormap=COLORMAP,color_limits=COLOR_LIMITS,drawer=draw,publish_options...)
    pair=read_pair(input_dir);s=read_summary(summary_path,input_dir)
    names=("temperature_comparison.png","boundary_heat.png","burgers_fields.png","convergence.png","plots.toml")
    staged(output_dir,names;publish_options...) do stage
        display=drawer(stage,pair,s;colormap,color_limits)
        doc=Dict("schema_version"=>1,"task_id"=>"N07","run_id"=>s["run_id"],"source_sha256"=>input_hashes(input_dir),"source_summary_sha256"=>file_sha(summary_path),"display_grid"=>DISPLAY_GRID,"display_time"=>1.,"cases"=>sort(collect(keys(s["cases"]))),"colormap"=>string(colormap),"temperature_color_limits"=>collect(color_limits),"error_colormap"=>"RdBu","burgers_error_limits"=>display.errors,"burgers_field_limits"=>display.field_limits,"residual_display_scale"=>1e-15,"figure_sha256"=>Dict(n=>file_sha(joinpath(stage,n)) for n in names if endswith(n,".png")))
        write_toml(joinpath(stage,"plots.toml"),doc)
    end
    println("N07を保存データから再作図しました。")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
