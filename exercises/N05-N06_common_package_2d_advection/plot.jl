module N06Plots
include("provided_support.jl")
using Plots
# 学生の編集対象: 物理条件・診断値を変えずに表示だけを調整する。
const COLORMAP = :viridis
const COLOR_LIMITS = (0.4,1.6)
const DISPLAY_CASE = "n080x060"
function draw(stage,data,s;colormap=COLORMAP,color_limits=COLOR_LIMITS,display_case=DISPLAY_CASE)
    require(length(color_limits)==2 && all(isfinite,color_limits) && color_limits[1]<color_limits[2],"色範囲が不正です")
    require(haskey(data.cases,display_case),"表示ケースがありません: $display_case")
    c=data.cases[display_case]; time=last(c["time"]); initial=c["u"][:,:,1]; final=c["u"][:,:,end]
    exact=exact_field(c["x"],c["y"],time;cx=data.metadata["cx"],cy=data.metadata["cy"]); error=final-exact
    panels=[]
    for (u,title) in ((initial,"Initial, t=0"),(final,"Numerical, t=$time"),(exact,"Exact, t=$time"),(error,"Numerical - exact"))
        iserror=u===error; amplitude=max(maximum(abs,error),eps())
        push!(panels,heatmap(c["x"],c["y"],permutedims(u);xlabel="x",ylabel="y",title,
            color=iserror ? :RdBu : colormap,clims=iserror ? (-amplitude,amplitude) : color_limits,
            xlims=(0.,2.),ylims=(0.,1.),aspect_ratio=:equal))
    end
    savefig(plot(panels...;layout=(2,2),size=(1100,660),margin=5Plots.mm),joinpath(stage,"fields.png"))
    mass=plot(;xlabel="Saved time",ylabel="Mass change / 1e-15",legend=:bottomleft)
    errors=plot(;xlabel="Saved time",ylabel="L2 error",legend=:topleft)
    for id in s["convergence"]["case_ids"]
        d=s["cases"][id]
        plot!(mass,d["time"],(d["mass"].-first(d["mass"]))./1e-15;label=id,marker=:circle)
        plot!(errors,d["time"],d["l2_error"];label=id,marker=:circle)
    end
    savefig(plot(mass,errors;layout=(1,2),size=(1000,400),margin=5Plots.mm),joinpath(stage,"diagnostics.png"))
    ids=s["convergence"]["case_ids"]; h=[s["cases"][id]["dx"] for id in ids]; err=s["convergence"]["errors"]
    p=plot(h,err;xscale=:log10,yscale=:log10,xlabel="dx (dy refined together)",ylabel="Final L2 error",label="Numerical",marker=:circle,legend=:topleft,size=(700,460))
    plot!(p,h,first(err).*h./first(h);label="First order",linestyle=:dash)
    savefig(p,joinpath(stage,"convergence.png"))
end
function main(;input_path=joinpath(DEFAULT_OUTPUT_DIR,FIELD_NAME),summary_path=joinpath(DEFAULT_OUTPUT_DIR,"summary.toml"),output_dir=DEFAULT_OUTPUT_DIR,
    colormap=COLORMAP,color_limits=COLOR_LIMITS,display_case=DISPLAY_CASE,publish_options...)
    data=read_fields(input_path); s=read_summary(summary_path,input_path)
    names=("fields.png","diagnostics.png","convergence.png","plots.toml")
    staged(output_dir,names;publish_options...) do stage
        draw(stage,data,s;colormap,color_limits,display_case)
        write_toml(joinpath(stage,"plots.toml"),Dict("schema_version"=>1,"source_fields_sha256"=>file_sha(input_path),
            "source_summary_sha256"=>file_sha(summary_path),"display_case"=>display_case,"time"=>[first(data.cases[display_case]["time"]),last(data.cases[display_case]["time"])],
            "colormap"=>string(colormap),"color_limits"=>collect(color_limits),"error_colormap"=>"RdBu","mass_display_scale"=>1e-15))
    end
    println("N06を保存場から再作図しました。")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
