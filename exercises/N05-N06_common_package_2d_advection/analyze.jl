module N06Analysis
include("provided_support.jl")
"""等幅周期格子の空間分散。未完成はnothing、完成後は有限・非負のFloat64。"""
function spatial_variance(u)
    # TODO(N06): 保存場の平均を使って空間分散を求める。
    return nothing
end
function main(;input_path=joinpath(DEFAULT_OUTPUT_DIR,FIELD_NAME),output_dir=DEFAULT_OUTPUT_DIR,variance=spatial_variance,publish_options...)
    data=read_fields(input_path)
    doc=diagnostics(data,variance); doc["source_fields_sha256"]=file_sha(input_path)
    staged(output_dir,("summary.toml",);publish_options...) do stage
        write_toml(joinpath(stage,"summary.toml"),doc)
    end
    println(doc["diagnostics_complete"] ? "N06の基本診断と分散を保存しました。plot.jlを再実行してください。" : "N06の基本診断を保存しました。spatial_varianceは未完成です（分散は未保存）。")
    doc
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
