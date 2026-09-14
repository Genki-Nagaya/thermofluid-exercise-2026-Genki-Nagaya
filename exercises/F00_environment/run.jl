if !isdefined(Main, :CourseWorkflow)
    include(joinpath(@__DIR__, "..", "..", "scripts", "lib", "CourseWorkflow.jl"))
end

module F00Environment

using Main.CourseWorkflow

export ObservedCheck,
    PreflightReport,
    REQUIRED_JULIA_VERSION,
    SUPPORTED_AGENTS,
    classify_runtime,
    collect_preflight,
    parse_preflight_arguments,
    print_preflight,
    run_f00_preflight

const REQUIRED_JULIA_VERSION = v"1.13.0"
const SUPPORTED_AGENTS = ("copilot", "codex", "amazon-q")

struct ObservedCheck
    id::Symbol
    passed::Bool
    observed::String
    action::String
end

struct PreflightReport
    runtime::ObservedCheck
    workspace::ObservedCheck
    julia::ObservedCheck
    git::ObservedCheck
end

function classify_runtime(platform, kernel_release)
    platform_name = Symbol(platform)
    release = strip(String(kernel_release))
    lower_release = lowercase(release)

    if platform_name == :windows
        return (
            kind = :native_windows,
            passed = false,
            observed = "native Windows",
            action = "WindowsではWSL2 Ubuntu 24.04を起動し、Linux側のJulia・Git・SSH・agentを使用してください。",
        )
    elseif platform_name == :macos
        return (
            kind = :macos,
            passed = true,
            observed = "native macOS",
            action = "",
        )
    elseif platform_name == :linux
        isempty(release) &&
            return (
                kind = :unknown_linux,
                passed = false,
                observed = "Linuxのkernel releaseを判定できません",
                action = "WSL2またはnative Linuxの端末でkernel情報を確認し、F00を再実行してください。",
            )
        if occursin("microsoft", lower_release) && occursin("wsl2", lower_release)
            return (
                kind = :wsl2,
                passed = true,
                observed = "WSL2 Ubuntu 24.04相当 (kernel: $release)",
                action = "",
            )
        elseif occursin("microsoft", lower_release)
            return (
                kind = :wsl1,
                passed = false,
                observed = "WSL1相当 (kernel: $release)",
                action = "WSL2 Ubuntu 24.04へ更新してから、Linux側でF00を再実行してください。",
            )
        else
            return (
                kind = :native_linux,
                passed = true,
                observed = "native Linux (kernel: $release)",
                action = "",
            )
        end
    end

    (
        kind = :unknown,
        passed = false,
        observed = "未対応の実行環境: $platform_name",
        action = "WSL2 Ubuntu 24.04、native macOS、またはnative LinuxでF00を再実行してください。",
    )
end

default_platform_probe() =
    Sys.iswindows() ? :windows :
    Sys.isapple() ? :macos :
    Sys.islinux() ? :linux :
    :unknown

function default_kernel_probe()
    if Sys.islinux()
        path = "/proc/sys/kernel/osrelease"
        return isfile(path) ? strip(read(path, String)) : ""
    end
    ""
end

default_runtime_probe() = (
    platform = default_platform_probe(),
    kernel_release = default_kernel_probe(),
    workspace = pwd(),
)

function workspace_check(runtime, workspace)
    path = normpath(String(workspace))
    linux_home_workspace =
        startswith(path, "/home/") && length(path) > length("/home/")
    outside_linux_home = runtime.kind == :wsl2 && !linux_home_workspace
    passed = runtime.passed && !outside_linux_home
    observed = "pwd: $path"

    if outside_linux_home
        if path == "/mnt/c" || startswith(path, "/mnt/c/")
            observed *= " (Windows側/mnt/c)"
        else
            observed *= " (WSL2のLinux filesystem外)"
        end
        action = "WSL2では学生リポジトリを /home/<user>/... にcloneし、Linux filesystem内でF00を実行してください。"
    elseif !runtime.passed
        action = runtime.action
    else
        action = ""
    end
    ObservedCheck(:workspace, passed, observed, action)
end

function default_command_probe(program, arguments)
    executable = Sys.which(program)
    isnothing(executable) &&
        return (available = false, detail = "$(program)がPATH上に見つかりません")

    stdout = IOBuffer()
    stderr = IOBuffer()
    command = Cmd([executable, arguments...])
    process = run(pipeline(ignorestatus(command), stdout = stdout, stderr = stderr))
    output =
        strip(join(filter(!isempty, [String(take!(stdout)), String(take!(stderr))]),
))
    detail = isempty(output) ? "$(program)の終了コード: $(process.exitcode)" : output
    (available = process.exitcode == 0, detail = detail)
end

function collect_preflight(;
    version_probe = () -> VERSION,
    command_probe = default_command_probe,
    runtime_probe = default_runtime_probe,
)
    runtime_input = runtime_probe()
    runtime_result = classify_runtime(runtime_input.platform, runtime_input.kernel_release)
    runtime_check = ObservedCheck(
        :runtime,
        runtime_result.passed,
        runtime_result.observed,
        runtime_result.action,
    )
    workspace = workspace_check(runtime_result, runtime_input.workspace)

    julia_version = version_probe()
    julia_check = ObservedCheck(
        :julia,
        julia_version == REQUIRED_JULIA_VERSION,
        string(julia_version),
        "JuliaupでJulia 1.13.0をインストールして選択し、この確認を再実行してください。",
    )

    git_probe = command_probe("git", ["--version"])
    git_check = ObservedCheck(
        :git,
        git_probe.available,
        String(git_probe.detail),
        "Gitをインストールし、GitコマンドをPATHから実行できることを確認してください。",
    )

    PreflightReport(runtime_check, workspace, julia_check, git_check)
end

function parse_preflight_arguments(arguments)
    vscode_confirmed = false
    github_confirmed = false
    agent = nothing
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        if argument == "--confirm-vscode"
            vscode_confirmed &&
                throw(ArgumentError("--confirm-vscodeは1回だけ指定できます"))
            vscode_confirmed = true
            index += 1
        elseif argument == "--confirm-github"
            github_confirmed &&
                throw(ArgumentError("--confirm-githubは1回だけ指定できます"))
            github_confirmed = true
            index += 1
        elseif argument == "--confirm-agent"
            isnothing(agent) || throw(ArgumentError("--confirm-agentは1回だけ指定できます"))
            index == length(arguments) &&
                throw(ArgumentError("--confirm-agentには製品名が必要です"))
            candidate = arguments[index + 1]
            candidate in SUPPORTED_AGENTS || throw(
                ArgumentError(
                    "未対応のAIエージェント「$candidate」です。copilot、codex、amazon-qから選んでください",
                ),
            )
            agent = candidate
            index += 2
        else
            throw(ArgumentError("「preflight」の不明な引数です: $argument"))
        end
    end
    (; vscode_confirmed, github_confirmed, agent)
end

function print_observed_check(io, label, check)
    status = check.passed ? "PASS" : "NEEDS SETUP"
    println(io, "  [$status] $label: $(check.observed)")
    check.passed || println(io, "    対応: $(check.action)")
end

function print_preflight(
    io,
    report;
    vscode_confirmed = false,
    github_confirmed = false,
    agent = nothing,
)
    println(io, "端末で確認した項目")
    print_observed_check(io, "Runtime", report.runtime)
    print_observed_check(io, "作業ディレクトリ (pwd)", report.workspace)
    print_observed_check(io, "Julia", report.julia)
    print_observed_check(io, "Git", report.git)
    println(io)
    println(io, "手動確認")
    println(
        io,
        "  [$(vscode_confirmed ? "CONFIRMED" : "NOT CONFIRMED")] VS Code・Julia拡張機能 (WindowsはRemote - WSL)",
    )
    println(
        io,
        "  [$(github_confirmed ? "CONFIRMED" : "NOT CONFIRMED")] GitHubへのサインインとリポジトリへのアクセス",
    )
    agent_label = isnothing(agent) ? "未指定" : agent
    println(
        io,
        "  [$(isnothing(agent) ? "NOT CONFIRMED" : "CONFIRMED")] 正式対応AIエージェント: $agent_label",
    )
    nothing
end

function run_f00_preflight(
    root;
    report = collect_preflight(),
    vscode_confirmed = false,
    github_confirmed = false,
    agent = nothing,
    persist_progress = save_progress,
    io = stdout,
)
    !isnothing(agent) &&
        !(agent in SUPPORTED_AGENTS) &&
        throw(ArgumentError("未対応のAIエージェント「$agent」です"))
    print_preflight(io, report; vscode_confirmed, github_confirmed, agent)

    observed_pass =
        report.runtime.passed &&
        report.workspace.passed &&
        report.julia.passed &&
        report.git.passed
    manual_pass = vscode_confirmed && github_confirmed && !isnothing(agent)
    if !(observed_pass && manual_pass)
        println(io)
        println(
            io,
            "F00の進捗は更新されませんでした。すべての機械観測と3項目の手動確認を完了してください。",
        )
        return false
    end

    progress_path = joinpath(root, "course_progress.toml")
    state = load_progress(progress_path)
    if state.current == "F01" && state.completed == ["F00"]
        println(io)
        println(io, "F00は完了済みです。現在の課題はF01のままです。")
        return true
    end
    state.current == "F00" && isempty(state.completed) ||
        throw(ArgumentError("F00の事前診断は初期F00進捗だけを更新できます"))

    advanced = ProgressState(state.schema_version, state.ordered, ["F00"], "F01")
    persist_progress(progress_path, advanced)
    println(io)
    println(
        io,
        "F00が完了しました。現在の課題はF01です。F00用のbranchやPRは作成しないでください。",
    )
    true
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    F00Environment.print_preflight(stdout, F00Environment.collect_preflight())
    println()
    println(
        "このスクリプトは診断専用です。https://t2lab-it.github.io/thermofluid-exercise-2026/assignments/F00.html に従い、scripts/course.jlからF00を完了してください。",
    )
end
