using Test
using ThermofluidExercise

isempty(ARGS) || error("使い方: julia --project=. -e 'using Pkg; Pkg.test()'")
const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(REPO_ROOT, "scripts", "lib", "CourseWorkflow.jl"))
include(joinpath(REPO_ROOT, "scripts", "lib", "ResultLimits.jl"))
include(joinpath(REPO_ROOT, "scripts", "lib", "CourseTests.jl"))
using .CourseTests

run_course_tests(REPO_ROOT) || error("現在課題または共通検査が失敗しました。上の失敗詳細を確認してください")
