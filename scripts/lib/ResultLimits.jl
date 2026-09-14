module ResultLimits

using ..CourseWorkflow: UNIT_DIRECTORIES, TASK_IDS_BY_UNIT
export FILE_LIMIT, TASK_LIMIT, TOTAL_LIMIT, check_result_limits

const FILE_LIMIT = 5 * 1024^2
const TASK_LIMIT = 10 * 1024^2
const TOTAL_LIMIT = 100 * 1024^2

function result_roots(root)
    roots = String[]
    exercises = joinpath(root, "exercises")
    if isdir(exercises)
        for (directory, children, _) in walkdir(exercises)
            if "results" in children
                push!(roots, joinpath(directory, "results"))
                filter!(!=("results"), children)
            end
        end
    end
    # An old output tree must not silently escape the size checks during migration.
    legacy = joinpath(root, "results")
    isdir(legacy) && push!(roots, legacy)
    roots
end

function content_id(root, path)
    parts = splitpath(relpath(path, root))
    length(parts) >= 4 && parts[1] == "exercises" && parts[3] == "results" || return nothing
    unit = findfirst(==(parts[2]), UNIT_DIRECTORIES)
    isnothing(unit) && return nothing
    ids = TASK_IDS_BY_UNIT[unit]
    length(ids) == 1 && return only(ids)
    length(parts) >= 5 && parts[4] in ids ? parts[4] : nothing
end

function check_result_limits(root; file_limit=FILE_LIMIT, task_limit=TASK_LIMIT, total_limit=TOTAL_LIMIT)
    for (name, limit) in (("file_limit", file_limit), ("task_limit", task_limit), ("total_limit", total_limit))
        limit >= 0 || throw(ArgumentError("`$name`は0以上にしてください"))
    end
    isdir(root) || throw(ArgumentError("リポジトリルートがディレクトリではありません: $root"))
    violations = String[]
    sizes = Dict{String,Int}()
    total = 0
    for results in result_roots(root), (directory, _, files) in walkdir(results), name in sort(files)
        path = joinpath(directory, name)
        relative = relpath(path, root)
        bytes = filesize(path)
        total += bytes
        id = content_id(root, path)
        if isnothing(id)
            push!(violations, "出力の配置を確認してください: $relative（$(bytes)バイト、ファイル上限$(file_limit)バイト）")
        else
            sizes[id] = get(sizes, id, 0) + bytes
        end
        bytes > file_limit && push!(violations,
            "ファイル上限を超えました: $relative は$(bytes)バイトです（上限$(file_limit)バイト）")
    end
    for id in sort!(collect(keys(sizes)))
        sizes[id] > task_limit && push!(violations,
            "課題上限を超えました: $id のresults/合計は$(sizes[id])バイトです（上限$(task_limit)バイト）")
    end
    total > total_limit && push!(violations,
        "成果物全体の上限を超えました: results/合計は$(total)バイトです（上限$(total_limit)バイト）")
    violations
end

end
