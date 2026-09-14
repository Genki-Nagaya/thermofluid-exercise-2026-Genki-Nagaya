using Test

const F00_REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const F00_RUN_SCRIPT = joinpath(F00_REPO_ROOT, "exercises", "F00_environment", "run.jl")

if !isdefined(Main, :F00Environment)
    include(F00_RUN_SCRIPT)
end
if !isdefined(Main, :CourseWorkflow)
    include(joinpath(F00_REPO_ROOT, "scripts", "lib", "CourseWorkflow.jl"))
end

using .F00Environment
using .CourseWorkflow

const WSL2_KERNEL = "5.15.167.4-microsoft-standard-WSL2"

contains_japanese(text::AbstractString) = occursin(r"[ぁ-んァ-ヶ一-龠]", text)

function runtime_input(runtime_kind, workspace)
    platform, kernel_release =
        if runtime_kind == :wsl2
            (:linux, WSL2_KERNEL)
        elseif runtime_kind == :native_linux
            (:linux, "6.8.0-generic")
        elseif runtime_kind == :macos
            (:macos, "Darwin Kernel Version 24.0.0")
        elseif runtime_kind == :windows
            (:windows, "")
        elseif runtime_kind == :wsl1
            (:linux, "4.4.0-Microsoft")
        elseif runtime_kind == :unknown_linux
            (:linux, "")
        else
            error("unknown runtime fixture: $runtime_kind")
        end
    (; platform, kernel_release, workspace)
end

function f00_report(;
    runtime_kind=:native_linux,
    workspace="/home/student/thermofluid-exercise",
    julia_ok=true,
    git_ok=true,
)
    versions = () -> julia_ok ? v"1.13.0" : v"1.12.6"
    commands = function (program, arguments)
        program == "git" ||
            error("unexpected command probe: $program $(join(arguments, ' '))")
        (available=git_ok, detail=git_ok ? "git version test" : "git not found")
    end
    runtime_probe = () -> runtime_input(runtime_kind, workspace)
    collect_preflight(; version_probe=versions, command_probe=commands, runtime_probe)
end

function initial_f00_repo()
    root = mktempdir()
    state = ProgressState(2, ORDERED_UNITS, String[], "F00")
    save_progress(joinpath(root, "course_progress.toml"), state)
    root
end

@testset "F00 environment preflight" begin
    @testset "runtime classification is injectable and explicit" begin
        @test classify_runtime(:linux, WSL2_KERNEL).kind == :wsl2
        @test classify_runtime(:linux, WSL2_KERNEL).passed
        @test classify_runtime(:linux, "6.8.0-generic").kind == :native_linux
        @test classify_runtime(:linux, "6.8.0-generic").passed
        @test classify_runtime(:macos, "Darwin").kind == :macos
        @test classify_runtime(:macos, "Darwin").passed
        @test classify_runtime(:windows, "").kind == :native_windows
        @test !classify_runtime(:windows, "").passed
        @test classify_runtime(:linux, "4.4.0-Microsoft").kind == :wsl1
        @test !classify_runtime(:linux, "4.4.0-Microsoft").passed
        @test classify_runtime(:linux, "").kind == :unknown_linux
        @test !classify_runtime(:linux, "").passed

        @test f00_report(runtime_kind=:wsl2).runtime.passed
        @test f00_report(runtime_kind=:macos).runtime.passed
        @test !f00_report(runtime_kind=:windows).runtime.passed
        @test !f00_report(runtime_kind=:wsl1).runtime.passed
        @test !f00_report(runtime_kind=:unknown_linux).runtime.passed

        mounted = f00_report(runtime_kind=:wsl2, workspace="/mnt/c/course")
        @test !mounted.workspace.passed
        @test occursin("/home/<user>", mounted.workspace.action)

        @test f00_report(runtime_kind=:wsl2, workspace="/home/student/course").workspace.passed
        for workspace in ("/mnt/c/course", "/mnt/d/course", "/tmp/course")
            outside_linux_home = f00_report(runtime_kind=:wsl2, workspace=workspace)
            @test !outside_linux_home.workspace.passed
            @test occursin("/home/<user>", outside_linux_home.workspace.action)
        end
    end

    @testset "machine-observed checks are structured and exact" begin
        report = f00_report()
        @test report.runtime.id == :runtime
        @test report.runtime.passed
        @test report.workspace.id == :workspace
        @test report.workspace.passed
        @test report.julia.id == :julia
        @test report.julia.passed
        @test report.julia.observed == "1.13.0"
        @test report.git.id == :git
        @test report.git.passed

        wrong_patch = f00_report(julia_ok=false)
        @test !wrong_patch.julia.passed
        @test contains_japanese(wrong_patch.julia.action)
        @test occursin("1.13.0", wrong_patch.julia.action)

        missing_git = f00_report(git_ok=false)
        @test !missing_git.git.passed
        @test contains_japanese(missing_git.git.action)
        @test occursin("Git", missing_git.git.action)

        missing_workspace = f00_report(runtime_kind=:wsl2, workspace="/mnt/c/course")
        @test !missing_workspace.workspace.passed
        @test contains_japanese(missing_workspace.workspace.action)
        @test occursin("/mnt/c", missing_workspace.workspace.observed)
    end

    @testset "manual confirmations are explicit and unambiguous" begin
        @test parse_preflight_arguments(String[]) ==
            (vscode_confirmed=false, github_confirmed=false, agent=nothing)
        @test parse_preflight_arguments([
            "--confirm-vscode", "--confirm-github", "--confirm-agent", "copilot",
        ]) == (vscode_confirmed=true, github_confirmed=true, agent="copilot")
        @test parse_preflight_arguments([
            "--confirm-agent", "codex", "--confirm-github", "--confirm-vscode",
        ]) == (vscode_confirmed=true, github_confirmed=true, agent="codex")
        @test parse_preflight_arguments([
            "--confirm-vscode", "--confirm-github", "--confirm-agent", "amazon-q",
        ]) == (vscode_confirmed=true, github_confirmed=true, agent="amazon-q")
        @test_throws ArgumentError parse_preflight_arguments(["--confirm-agent", "claude"])
        @test_throws ArgumentError parse_preflight_arguments([
            "--confirm-agent", "copilot", "--confirm-agent", "codex",
        ])
        @test_throws ArgumentError parse_preflight_arguments(["--confirm-github", "--confirm-github"])
        @test_throws ArgumentError parse_preflight_arguments(["--confirm-vscode", "--confirm-vscode"])
        @test_throws ArgumentError parse_preflight_arguments(["--unknown"])
    end

    @testset "progress changes only after observed and manual gates pass" begin
        root = initial_f00_repo()
        progress_path = joinpath(root, "course_progress.toml")
        before = read(progress_path, String)

        output = IOBuffer()
        @test !run_f00_preflight(root; report=f00_report(), io=output)
        @test read(progress_path, String) == before
        text = String(take!(output))
        @test contains_japanese(text)
        for identifier in ("Runtime", "pwd", "Julia", "Git", "VS Code", "F00")
            @test occursin(identifier, text)
        end

        output = IOBuffer()
        @test !run_f00_preflight(
            root;
            report=f00_report(),
            github_confirmed=true,
            io=output,
        )
        @test read(progress_path, String) == before

        for failed_report in (
            f00_report(julia_ok=false),
            f00_report(git_ok=false),
            f00_report(runtime_kind=:windows),
            f00_report(runtime_kind=:wsl2, workspace="/mnt/c/course"),
        )
            output = IOBuffer()
            @test !run_f00_preflight(
                root;
                report=failed_report,
                vscode_confirmed=true,
                github_confirmed=true,
                agent="codex",
                io=output,
            )
            @test read(progress_path, String) == before
        end

        output = IOBuffer()
        @test run_f00_preflight(
            root;
            report=f00_report(),
            vscode_confirmed=true,
            github_confirmed=true,
            agent="codex",
            io=output,
        )
        state = load_progress(progress_path)
        @test state.completed == ["F00"]
        @test state.current == "F01"
        completed_output = String(take!(output))
        @test contains_japanese(completed_output)
        @test occursin("F00", completed_output)
        @test occursin("F01", completed_output)
    end

    @testset "completion is atomic and idempotent without Git side effects" begin
        root = initial_f00_repo()
        progress_path = joinpath(root, "course_progress.toml")
        before = read(progress_path, String)
        writes = Ref(0)
        failing_persistence = function (path, state)
            writes[] += 1
            error("injected persistence failure")
        end
        @test_throws ErrorException run_f00_preflight(
            root;
            report=f00_report(),
            vscode_confirmed=true,
            github_confirmed=true,
            agent="copilot",
            persist_progress=failing_persistence,
            io=IOBuffer(),
        )
        @test writes[] == 1
        @test read(progress_path, String) == before

        @test run_f00_preflight(
            root;
            report=f00_report(),
            vscode_confirmed=true,
            github_confirmed=true,
            agent="copilot",
            io=IOBuffer(),
        )
        after = read(progress_path, String)
        writes[] = 0
        @test run_f00_preflight(
            root;
            report=f00_report(),
            vscode_confirmed=true,
            github_confirmed=true,
            agent="amazon-q",
            persist_progress=(path, state) -> (writes[] += 1),
            io=IOBuffer(),
        )
        @test writes[] == 0
        @test read(progress_path, String) == after
        @test !isdir(joinpath(root, "learning_logs"))
        @test !isdir(joinpath(root, ".git"))
    end
end

function f00_command_result(command)
    stdout = IOBuffer()
    stderr = IOBuffer()
    process = run(pipeline(ignorestatus(command), stdout=stdout, stderr=stderr))
    (
        exitcode=process.exitcode,
        stdout=String(take!(stdout)),
        stderr=String(take!(stderr)),
    )
end

function f00_cli_root()
    root = mktempdir()
    CourseWorkflow.save_progress(
        joinpath(root, "course_progress.toml"),
        CourseWorkflow.ProgressState(2, CourseWorkflow.ORDERED_UNITS, String[], "F00"),
    )
    root
end

@testset "F00 course CLI wiring" begin
    help_command = Cmd(Cmd([
        Base.julia_cmd().exec...,
        "--startup-file=no",
        "--project=$F00_REPO_ROOT",
        joinpath(F00_REPO_ROOT, "scripts", "course.jl"),
        "--help",
    ]); dir=F00_REPO_ROOT)
    help_result = f00_command_result(help_command)
    @test help_result.exitcode == 0
    @test occursin("--confirm-vscode", help_result.stdout)
    @test occursin("--confirm-github", help_result.stdout)
    @test occursin("--confirm-agent", help_result.stdout)

    process_root = f00_cli_root()
    progress_path = joinpath(process_root, "course_progress.toml")
    before = read(progress_path, String)
    empty_path = mktempdir()
    command = Cmd(Cmd([
        Base.julia_cmd().exec...,
        "--startup-file=no",
        "--project=$F00_REPO_ROOT",
        joinpath(F00_REPO_ROOT, "scripts", "course.jl"),
        "preflight",
        "--confirm-vscode",
        "--confirm-github",
        "--confirm-agent",
        "codex",
    ]); dir=process_root)
    result = f00_command_result(addenv(command, "PATH" => empty_path))
    @test result.exitcode != 0
    @test contains_japanese(result.stdout)
    @test occursin("NEEDS SETUP", result.stdout)
    @test occursin("F00", result.stdout)
    @test read(progress_path, String) == before
end
