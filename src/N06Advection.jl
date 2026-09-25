module N06
import ..ThermofluidExercise as Common
export stable_timestep, advection_step!
# 提供: 入力検証。数値式は下のTODO 2か所で実装する。
function positive(v,name)
    v isa Real && isfinite(v) && v>0 || throw(ArgumentError("$name は有限な正値です"))
end
function speeds(cx,cy)
    all(v->v isa Real && isfinite(v) && v>=0,(cx,cy)) || throw(ArgumentError("cx,cyは有限な非負値です"))
end
"""合成CFL上限。片方向ゼロ可、両方向ゼロ不可。"""
function stable_timestep(cx,cy,dx,dy;safety=0.8)
    speeds(cx,cy); positive(dx,"dx"); positive(dy,"dy"); positive(safety,"safety")
    safety<=1 && (cx>0 || cy>0) || throw(ArgumentError("safety<=1、少なくとも一方の速度は正です"))
    # TODO(N06): 方向別寄与を足した安定条件から刻みを求める。
    dt = error("未実装 N06: 合成時間刻み")
    positive(dt,"合成dt（表現範囲）")
    result=Float64(dt); positive(result,"Float64の合成dt")
    return result
end
"""u[i,j]はx[i],y[j]。周期全点を同じ旧行列から更新し、新旧は非alias。"""
function advection_step!(u_new,u_old,dt,dx,dy;cx=1.0,cy=0.5)
    u_new isa AbstractMatrix && u_old isa AbstractMatrix || throw(ArgumentError("新旧は行列です"))
    Common.validate_buffers(u_new,u_old)
    speeds(cx,cy); positive(dt,"dt"); positive(dx,"dx"); positive(dy,"dy")
    Cx=cx*(dt/dx); Cy=cy*(dt/dy)
    all(isfinite,(Cx,Cy)) && Cx+Cy<=1+32eps(Float64) || throw(ArgumentError("合成CFL<=1を超えています: Cx=$Cx, Cy=$Cy"))
    nx,ny=size(u_old)
    # TODO(N06): 両方向の周期隣接と流束を使い、全点を旧配列だけから更新する。
    error("未実装 N06: 周期二次元更新")
    return u_new
end
end
