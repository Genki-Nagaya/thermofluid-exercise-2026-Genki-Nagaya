module F01FirstPullRequest

export main, student_greeting

"""
    student_greeting(name::AbstractString) -> String

`name`の前後の空白を除き、`"Hello, <name>!"`を返す。
空の名前は無効とする。F01では印を付けたreturn文を完成させる。
"""
function student_greeting(name::AbstractString)::String
    normalized_name = strip(name)
    isempty(normalized_name) && throw(ArgumentError("名前を空にはできません"))

    # TODO(F01): この仮実装を`Hello, <normalized_name>!`へ置き換える。
    "TODO: implement student_greeting"
end

function main(name::AbstractString = "student"; io = stdout)
    message = student_greeting(name)
    println(io, message)
    message
end

end


if abspath(PROGRAM_FILE) == @__FILE__
    name = isempty(ARGS) ? "student" : only(ARGS)
    F01FirstPullRequest.main(name)
end
