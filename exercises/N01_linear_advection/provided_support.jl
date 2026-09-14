# N01の実行時入力検証と公式出力の詳細を担当する教材提供ファイルです。
# 数値計算の学習対象ではなく、受講生が読解・編集する必要はありません。

using Plots
using TOML

const DEFAULT_OUTPUT_DIR = joinpath(@__DIR__, "results")

"""初期条件の引数がN01の実行時の入力条件を満たすか確認する。"""
function validate_initial_condition_inputs(x, base, plateau, plateau_start, plateau_end)
    isempty(x) && throw(ArgumentError("xを空にすることはできません"))
    all(isfinite, x) || throw(ArgumentError("xの全要素を有限値にしてください"))
    all(isfinite, (base, plateau, plateau_start, plateau_end)) ||
        throw(ArgumentError("初期条件のパラメータを有限値にしてください"))
    all(diff(x) .> 0) || throw(ArgumentError("xを狭義単調増加にしてください"))
    first(x) <= plateau_start <= plateau_end <= last(x) ||
        throw(ArgumentError("矩形領域を計算領域内に置いてください"))
    return nothing
end

"""1ステップ更新の引数がN01の実行時の入力条件を満たすか確認する。"""
function validate_step_inputs(u_new, u_old, c, dt, dx)
    u_new === u_old && throw(ArgumentError("新旧で別々のバッファを使ってください"))
    length(u_new) == length(u_old) >= 3 ||
        throw(ArgumentError("二つのバッファを同じ長さの3点以上にしてください"))
    all(isfinite, u_old) || throw(ArgumentError("u_oldの全要素を有限値にしてください"))
    all(isfinite, (c, dt, dx)) || throw(ArgumentError("c、dt、dxを有限値にしてください"))
    c > 0 || throw(ArgumentError("N01では正の移流速度だけを扱います"))
    dt > 0 || throw(ArgumentError("dtは正にしてください"))
    dx > 0 || throw(ArgumentError("dxは正にしてください"))
    return nothing
end

"""境界条件の引数がN01の実行時の入力条件を満たすか確認する。"""
function validate_boundary_inputs(u, left_value)
    length(u) >= 2 || throw(ArgumentError("uには2点以上が必要です"))
    isfinite(left_value) || throw(ArgumentError("left_valueを有限値にしてください"))
    return nothing
end

"""時間発展の引数がN01の実行時の入力条件を満たすか確認する。"""
function validate_simulation_inputs(scheme, nx, c, cfl, t_final)
    scheme in (:upwind, :centered) ||
        throw(ArgumentError("schemeには:upwindまたは:centeredを指定してください"))
    nx isa Bool && throw(ArgumentError("nxには格子点数を表す整数を指定してください"))
    nx >= 3 || throw(ArgumentError("nxは3以上にしてください"))
    all(isfinite, (c, cfl, t_final)) ||
        throw(ArgumentError("c、cfl、t_finalを有限値にしてください"))
    c > 0 || throw(ArgumentError("N01では正の移流速度だけを扱います"))
    0 < cfl <= 1 || throw(ArgumentError("cflは0 < cfl <= 1を満たす必要があります"))
    t_final > 0 || throw(ArgumentError("t_finalは正にしてください"))
    return nothing
end

"""機械可読な診断量を書き出し、summary.tomlのパスを返す。"""
function write_summary(output_dir::AbstractString, upwind, centered)
    mkpath(output_dir)
    path = joinpath(output_dir, "summary.toml")
    summary = Dict(
        "course_id" => "N01",
        "grid" => Dict("nx" => length(upwind.x), "dx" => upwind.dx),
        "upwind" => summary_section("upwind-euler", upwind),
        "centered_euler" => summary_section("centered-euler", centered),
    )
    open(path, "w") do io
        TOML.print(io, summary; sorted = true)
    end
    return path
end

"""公式の比較図を二つ作り、それぞれのパスを返す。"""
function make_plots(output_dir::AbstractString, upwind, centered)
    mkpath(output_dir)
    upwind_path = joinpath(output_dir, "upwind.png")
    centered_path = joinpath(output_dir, "centered-euler.png")

    for (result, title, path) in (
        (upwind, "Upwind + Euler (stable)", upwind_path),
        (centered, "Centered + Euler (intentionally unstable)", centered_path),
    )
        plot(
            result.x,
            result.u0;
            label = "Initial",
            linewidth = 2,
            xlabel = "x",
            ylabel = "u",
            title = title,
        )
        plot!(result.x, result.u; label = "Final", linewidth = 2)
        savefig(path)
    end
    return (upwind = upwind_path, centered = centered_path)
end
