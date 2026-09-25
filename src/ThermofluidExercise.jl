module ThermofluidExercise

"""周期左隣。整数1 <= i <= n、Bool不可。"""
function periodic_left_index(i,n)
    error("未実装 N05: periodic_left_index")
end

"""周期右隣。整数1 <= i <= n、Bool不可。"""
function periodic_right_index(i,n)
    error("未実装 N05: periodic_right_index")
end

"""有限実数の線形流束。"""
function linear_flux(u,speed)
    error("未実装 N05: linear_flux")
end

"""有限実数のBurgers流束。負値も受け付ける。"""
function burgers_flux(u)
    error("未実装 N05: burgers_flux")
end

"""1始まり・浮動小数・同形状・各軸3点以上・非alias・有限旧値。"""
function validate_buffers(u_new,u_old)
    error("未実装 N05: validate_buffers")
end

"""正の有限値から (;steps,dt) を返し最終時刻へ合わせる。"""
function fit_timestep(dt_max,t_final)
    error("未実装 N05: fit_timestep")
end

module N01
"""N01の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function upwind_step!(args...;kwargs...)
    error("未実装 N05: N01.upwind_step!")
end
"""N01の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function centered_step!(args...;kwargs...)
    error("未実装 N05: N01.centered_step!")
end
"""N01の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function apply_boundary!(args...;kwargs...)
    error("未実装 N05: N01.apply_boundary!")
end
end
module N02
"""N02の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function nonlinear_upwind_step!(args...;kwargs...)
    error("未実装 N05: N02.nonlinear_upwind_step!")
end
"""N02の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function apply_boundary!(args...;kwargs...)
    error("未実装 N05: N02.apply_boundary!")
end
end
module N03
"""N03の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function diffusion_step!(args...;kwargs...)
    error("未実装 N05: N03.diffusion_step!")
end
"""N03の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function apply_boundary!(args...;kwargs...)
    error("未実装 N05: N03.apply_boundary!")
end
end
module N04
"""N04の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function stable_timestep(args...;kwargs...)
    error("未実装 N05: N04.stable_timestep")
end
"""N04の既存入口と同じ契約を保つ。N04更新ではmodelを明示する。"""
function advection_diffusion_step!(args...;kwargs...)
    error("未実装 N05: N04.advection_diffusion_step!")
end
end
include("N06Advection.jl")
include("N07Transport.jl")
end
