module N0506Run
include("N05.jl")
include("simulate.jl")
include("analyze.jl")
include("plot.jl")
function main(;output_dir=joinpath(@__DIR__,"results"),baseline_path=joinpath(output_dir,"N05","baseline.toml"),simulation=N06Simulation.main,analysis=N06Analysis.main,plotting=N06Plots.main,publish_options...)
    mktempdir() do stage
        N05Regression.verify(;baseline_path,output_dir=joinpath(stage,"N05"))
        target=joinpath(stage,"N06")
        simulation(;output_dir=target)
        analysis(;input_path=joinpath(target,"fields.h5"),output_dir=target)
        plotting(;input_path=joinpath(target,"fields.h5"),summary_path=joinpath(target,"summary.toml"),output_dir=target)
        names=vcat([joinpath("N05","regression.toml")],[joinpath("N06",n) for n in N06Simulation.OUTPUT_NAMES])
        N06Simulation.publish(stage,output_dir,names;publish_options...)
    end
    println("N05回帰とN06の全出力を反映しました。baselineは保持しています。")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
