module N07
import ..ThermofluidExercise as Common
export thermal_stable_timestep, thermal_fluxes!, thermal_step!, burgers_stable_timestep, burgers_step!
# 提供: 形状・境界・書込み前検証。数値処理はTODOで実装する。
require(ok,msg)=ok ? nothing : throw(ArgumentError(msg))
finite(v)=v isa Real && isfinite(v)
positive(v,name)=require(finite(v) && v>0,"$name は有限正値です")
const SIDES=(:west,:east,:south,:north)
function boundaries(kind::Symbol;value=0.0)
    b=(;kind,value=Float64(value)); bc=(;west=b,east=b,south=b,north=b)
    validate_boundary(bc,0.,0.); bc
end
function channel_boundaries(wall::Symbol)
    require(wall in (:dirichlet,:insulated),"上下壁は固定温度または断熱です")
    (;west=(;kind=:inflow,value=1.),east=(;kind=:outflow,value=0.),south=(;kind=wall,value=0.),north=(;kind=wall,value=0.))
end
function validate_boundary(bc,cx,cy)
    require(bc isa NamedTuple && keys(bc)==SIDES,"bcはwest,east,south,northです")
    for b in values(bc)
        require(b isa NamedTuple && keys(b)==(:kind,:value) && b.kind in (:periodic,:dirichlet,:insulated,:inflow,:outflow) && b.value isa Float64 && isfinite(b.value),"境界の種類または値が不正です")
        require(b.kind in (:dirichlet,:inflow) || b.value==0.,"値が不要な境界はvalue=0.0です")
    end
    for (a,b) in ((bc.west,bc.east),(bc.south,bc.north))
        require((a.kind==:periodic)==(b.kind==:periodic),"周期辺は対向辺で指定します")
    end
    kinds=map(b->b.kind,values(bc))
    if :inflow in kinds || :outflow in kinds
        require(kinds[1:2]==(:inflow,:outflow) && all(k->k in (:dirichlet,:insulated),kinds[3:4]) && cx>0 && cy==0,"開境界は左流入・右流出、cx>0,cy=0です")
    elseif !all(==(:periodic),kinds)
        require(all(k->k in (:periodic,:dirichlet,:insulated),kinds) && cx==cy==0,"閉じた非周期系の速度は0です")
    end
    nothing
end
function thermal_parameters(cx,cy,kappa,dx,dy,bc;safety=1.)
    require(all(finite,(cx,cy,kappa)) && kappa>=0,"速度は有限、kappaは有限非負です")
    positive(dx,"dx"); positive(dy,"dy"); positive(safety,"safety"); require(safety<=1,"safety<=1です")
    validate_boundary(bc,cx,cy)
end
function matrix(a;old=false)
    require(a isa AbstractMatrix{<:AbstractFloat} && all(>=(3),size(a)),"1始まり・各軸3以上の浮動小数行列です")
    require(axes(a)==map(Base.OneTo,size(a)),"1始まりの配列です")
    old && require(all(isfinite,a),"旧場は有限です")
end
function independent(arrays)
    for i in eachindex(arrays), j in i+1:length(arrays)
        require(!Base.mightalias(arrays[i],arrays[j]),"配列は互いに非aliasです")
    end
end
function thermal_buffers(new,old)
    matrix(new); matrix(old;old=true); require(size(new)==size(old),"新旧の形状が異なります"); independent((new,old))
end
function flux_buffers(T)
    matrix(T); nx,ny=size(T)
    (;adv_x=zeros(nx+1,ny),diff_x=zeros(nx+1,ny),adv_y=zeros(nx,ny+1),diff_y=zeros(nx,ny+1))
end
function validate_fluxes(f,T)
    matrix(T;old=true); nx,ny=size(T)
    require(f isa NamedTuple && keys(f)==(:adv_x,:diff_x,:adv_y,:diff_y),"流束配列名が不正です")
    for (a,shape) in zip(values(f),((nx+1,ny),(nx+1,ny),(nx,ny+1),(nx,ny+1)))
        require(a isa AbstractMatrix{<:AbstractFloat} && size(a)==shape && axes(a)==map(Base.OneTo,shape),"流束配列の形状・軸が不正です")
    end
    independent((T,values(f)...))
end
# 検証用の中心係数。学生はstable_timestepでこの上限の意味を説明して式を実装する。
function thermal_rate(cx,cy,kappa,dx,dy,bc)
    ax=any(b->b.kind in (:dirichlet,:inflow),(bc.west,bc.east)) ? 3 : 2
    ay=any(b->b.kind in (:dirichlet,:inflow),(bc.south,bc.north)) ? 3 : 2
    r=abs(cx)/dx+abs(cy)/dy+kappa*(ax/dx^2+ay/dy^2)
    require(isfinite(r),"安定条件が表現範囲外です"); r
end
function checked_timestep(dt)
    result=Float64(dt); positive(result,"dt"); result
end
function thermal_stable_timestep(cx,cy,kappa,dx,dy,bc;safety=.8)
    thermal_parameters(cx,cy,kappa,dx,dy,bc;safety)
    require(thermal_rate(cx,cy,kappa,dx,dy,bc)>0,"全係数0では刻みを決められません")
    # TODO(N07): 半セル境界を含む合成上限。
    dt=error("未実装 N07: 温度の安定刻み")
    checked_timestep(dt)
end
function thermal_fluxes!(fluxes,T,dx,dy;cx,cy,kappa,bc)
    thermal_parameters(cx,cy,kappa,dx,dy,bc); validate_fluxes(fluxes,T)
    # TODO(N07): 正x・正y方向、移流と拡散を分けて全ての面へ書く。
    error("未実装 N07: 温度の面流束")
    fluxes
end
function thermal_step!(Tnew,Told,dt,dx,dy;cx,cy,kappa,bc)
    thermal_parameters(cx,cy,kappa,dx,dy,bc); thermal_buffers(Tnew,Told); positive(dt,"dt")
    rate=thermal_rate(cx,cy,kappa,dx,dy,bc)
    require(dt*rate<=1+32eps(Float64),"温度の安定上限を超えています")
    # TODO(N07): Common.validate_buffersも再利用し、同じ旧場の面流束を一度構築して更新と辺別レートに使う。
    error("未実装 N07: 温度更新と境界熱流束")
end
function burgers_parameters(u,v,nu,dx,dy;safety=1.)
    matrix(u;old=true); matrix(v;old=true); require(size(u)==size(v),"二成分の形状が異なります"); independent((u,v))
    require(finite(nu) && nu>=0,"nuは有限非負です"); positive(dx,"dx"); positive(dy,"dy"); positive(safety,"safety"); require(safety<=1,"safety<=1です")
end
function burgers_rate(u,v,nu,dx,dy)
    r=maximum(abs(u[k])/dx+abs(v[k])/dy for k in eachindex(u,v))+2nu*(1/dx^2+1/dy^2)
    require(isfinite(r),"安定条件が表現範囲外です"); r
end
function burgers_stable_timestep(u,v,nu,dx,dy;safety=.8)
    burgers_parameters(u,v,nu,dx,dy;safety)
    require(burgers_rate(u,v,nu,dx,dy)>0,"全速度・粘性0では刻みを決められません")
    # TODO(N07): このステップの旧場の局所速度から合成刻みを決める。
    dt=error("未実装 N07: Burgersの動的刻み")
    checked_timestep(dt)
end
function burgers_step!(unew,vnew,uold,vold,dt,dx,dy;nu)
    burgers_parameters(uold,vold,nu,dx,dy); thermal_buffers(unew,uold); thermal_buffers(vnew,vold)
    independent((unew,vnew,uold,vold)); positive(dt,"dt")
    require(dt*burgers_rate(uold,vold,nu,dx,dy)<=1+32eps(Float64),"Burgersの安定上限を超えています")
    # TODO(N07): Common.validate_buffersと周期添字を再利用し、正負風上と拡散で二成分を同じ旧場から更新。
    error("未実装 N07: 二成分Burgers更新")
    (;u=unew,v=vnew)
end
end
