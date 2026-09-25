module N07Analysis
include("provided_support.jl")
"""辺別の区間積分から熱量変化を再構成する。流入を正とする。"""
function heat_budget(time,heat,advective_integrals,diffusive_integrals)
    validate_budget(time,heat,advective_integrals,diffusive_integrals)
    # TODO(N07): 各区間の正味流入、累積流入、初期熱量からの変化との差。
    error("未実装 N07: 熱収支解析")
end
function main(;input_dir=DEFAULT_OUTPUT_DIR,output_dir=DEFAULT_OUTPUT_DIR,budget=heat_budget,publish_options...)
    data=read_pair(input_dir);doc=diagnostics(data,budget);doc["source_sha256"]=input_hashes(input_dir)
    staged(output_dir,("summary.toml",);publish_options...) do stage
        write_toml(joinpath(stage,"summary.toml"),doc)
    end
    println("N07の保存場と区間熱輸送を解析しました。plot.jlを再実行してください。")
    doc
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
