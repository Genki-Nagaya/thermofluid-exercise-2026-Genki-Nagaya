module N07Run
include("simulate.jl")
include("analyze.jl")
include("plot.jl")
function main(;output_dir=N07Simulation.DEFAULT_OUTPUT_DIR,simulation=N07Simulation.main,analysis=N07Analysis.main,plotting=N07Plots.main,publish_options...)
    mktempdir() do stage
        simulation(;output_dir=stage)
        analysis(;input_dir=stage,output_dir=stage)
        plotting(;input_dir=stage,output_dir=stage)
        N07Simulation.check_complete(stage)
        N07Simulation.publish(stage,output_dir,N07Simulation.OUTPUT_NAMES;publish_options...)
    end
    println("N07の公式8出力をまとめて反映しました。")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
