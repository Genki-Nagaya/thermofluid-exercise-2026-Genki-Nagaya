using Test
using ThermofluidExercise

isempty(ARGS) || error("使い方: julia --project=. -e 'using Pkg; Pkg.test()'")
const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "scripts", "lib", "CourseWorkflow.jl"))
include(joinpath(REPO_ROOT, "scripts", "lib", "ResultLimits.jl"))
using .CourseWorkflow
using .ResultLimits
include(joinpath(REPO_ROOT, "test", "f00_preflight_test.jl"))

state = load_progress(joinpath(REPO_ROOT, "course_progress.toml"))
for unit in units_to_test(state)
    path = joinpath(REPO_ROOT, unit_directory(unit), "tests.jl")
    isfile(path) || error("$unit の必須テストがありません: $path")
    @testset "$unit / tests.jl" begin
        include(path)
    end
end
violations = check_result_limits(REPO_ROOT)
isempty(violations) || error(join(violations, '\n'))
