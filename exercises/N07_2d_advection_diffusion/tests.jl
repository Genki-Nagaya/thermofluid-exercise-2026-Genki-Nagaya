module N07Tests
using Test, ThermofluidExercise
const N=ThermofluidExercise.N07
include("provided_support.jl")
include("analyze.jl")
# 配布済みの必須テスト。自作欄と区別し、入力・期待値を保持する。
@testset "N07 independent faces, half-cell walls and heat" begin
    a=Float64[1 2 4 3;5 7 6 8;9 10 12 11]; saved=copy(a); b=similar(a)
    bc=N.channel_boundaries(:dirichlet); f=N.flux_buffers(a)
    @test N.thermal_stable_timestep(1.,0.,.05,.5,.25,bc)≈.8/5.0
    @test N.thermal_fluxes!(f,a,.5,.25;cx=1.,cy=0.,kappa=.05,bc)===f
    @test f.adv_x==[1 1 1 1;1 2 4 3;5 7 6 8;9 10 12 11]
    @test f.diff_x≈[0 -.2 -.6 -.4;-.4 -.5 -.2 -.5;-.4 -.3 -.6 -.3;0 0 0 0]
    @test f.adv_y==zeros(3,5)
    @test f.diff_y≈[-.4 -.2 -.4 .2 1.2;-2 -.4 .2 -.4 3.2;-3.6 -.2 -.4 .2 4.4]
    r=N.thermal_step!(b,a,.01,.5,.25;cx=1.,cy=0.,kappa=.05,bc)
    @test r.temperature===b
    @test collect(r.boundary_rates.advective)==[1.,-10.5,0.,0.]
    @test collect(r.boundary_rates.diffusive)≈[-.3,0.,-3.,-4.4]
    # Independent conservative sum using explicitly checked face arrays.
    want=a-.01*((f.adv_x[2:end,:]+f.diff_x[2:end,:]-f.adv_x[1:end-1,:]-f.diff_x[1:end-1,:])/.5+(f.diff_y[:,2:end]-f.diff_y[:,1:end-1])/.25)
    @test b≈want atol=1e-14
    @test .5*.25*sum(b-a)≈.01*(-17.2) atol=1e-14
    @test a==saved
    for wall in (:dirichlet,:insulated)
        bc=N.channel_boundaries(wall); q0=.125sum(a); total=0.; u=copy(a); v=similar(u)
        for _ in 1:20
            r=N.thermal_step!(v,u,.01,.5,.25;cx=1.,cy=0.,kappa=.05,bc)
            rates=r.boundary_rates
            @test .125sum(v-u)≈.01*(sum(rates.advective)+sum(rates.diffusive)) atol=1e-12
            wall==:insulated && @test rates.diffusive.south==rates.diffusive.north==0.
            total+=.01*(sum(rates.advective)+sum(rates.diffusive)); u,v=v,u
        end
        @test .125sum(u)-q0≈total atol=1e-12
    end
    for kind in (:periodic,:insulated,:dirichlet)
        bc=N.boundaries(kind;value=kind==:dirichlet ? 2. : 0.)
        N.thermal_step!(b,fill(2.,3,4),.01,.5,.25;cx=0.,cy=0.,kappa=.05,bc)
        @test b≈fill(2.,3,4)
    end
    bc=N.boundaries(:dirichlet); bound=1/(.05*(3/.5^2+3/.25^2))
    impulse=zeros(3,4); impulse[1,1]=1.
    N.thermal_step!(b,impulse,bound,.5,.25;cx=0.,cy=0.,kappa=.05,bc)
    @test minimum(b)>=-1e-15
    fill!(b,-99.)
    @test_throws ArgumentError N.thermal_step!(b,a,bound*(1+1e-8),.5,.25;cx=0.,cy=0.,kappa=.05,bc)
    @test all(==(-99.),b)
    @test_throws ArgumentError N.thermal_step!(b,a,1/(.1*(1/.5^2+1/.25^2)),.5,.25;cx=0.,cy=0.,kappa=.05,bc)
    # Signed periodic flow at both ends, old-state use and heat cancellation.
    bc=N.boundaries(:periodic)
    for (cx,cy,kappa) in ((-1.,.5,0.),(0.,0.,.05),(0.,0.,0.))
        r=N.thermal_step!(b,a,.01,.5,.25;cx,cy,kappa,bc)
        @test sum(b)≈sum(a) atol=1e-12
        @test sum(r.boundary_rates.advective)+sum(r.boundary_rates.diffusive)≈0. atol=1e-13
        if kappa==0 && cx!=0
            @test b[3,1]≈9-.01*((-1*(1-9))/.5+.5*(9-11)/.25)
        elseif cx==cy==kappa==0
            @test b==a
        end
    end
end
@testset "N07 two-component signed advective Burgers" begin
    u=[1. -2. 0. .5;-.5 2. -1. 0.;3. -1. .2 -2.]
    v=[-.5 1. 2. -1.;2. -1. 0. .5;0. .5 -2. 1.]
    originals=(copy(u),copy(v)); un=similar(u); vn=similar(v); dt=.001; dx=.5; dy=.25; nu=.05
    # Independent whole-array oracle uses circular translations, not production neighbor helpers.
    function expected(a)
        left=circshift(a,(1,0)); right=circshift(a,(-1,0)); down=circshift(a,(0,1)); up=circshift(a,(0,-1))
        ux=ifelse.(u.>=0,(a-left)/dx,(right-a)/dx)
        uy=ifelse.(v.>=0,(a-down)/dy,(up-a)/dy)
        a-dt*(u.*ux+v.*uy)+dt*nu*((left-2a+right)/dx^2+(down-2a+up)/dy^2)
    end
    r=N.burgers_step!(un,vn,u,v,dt,dx,dy;nu)
    @test r.u===un && r.v===vn
    @test un≈expected(u) atol=1e-14
    @test vn≈expected(v) atol=1e-14
    @test u==originals[1] && v==originals[2]
    @test N.burgers_stable_timestep(u,v,nu,dx,dy)≈.8/(9.0+2.0)
    @test N.burgers_stable_timestep(2u,2v,nu,dx,dy)<N.burgers_stable_timestep(u,v,nu,dx,dy)
    for a in (zeros(3,4),fill(2.,3,4)), b in (zeros(3,4),fill(-1.,3,4))
        N.burgers_step!(un,vn,a,b,.001,dx,dy;nu=0.)
        @test un==a && vn==b
    end
    N.burgers_step!(un,vn,u,zeros(3,4),dt,dx,dy;nu=0.); @test all(iszero,vn)
    arrays=[copy(u) for _ in 1:4]
    for i in 1:4,j in i+1:4
        args=copy(arrays); args[j]=args[i]
        @test_throws ArgumentError N.burgers_step!(args...,dt,dx,dy;nu)
    end
    fill!(un,-99.); fill!(vn,-88.)
    @test_throws ArgumentError N.burgers_step!(un,vn,u,v,1.,dx,dy;nu)
    @test all(==(-99.),un) && all(==(-88.),vn)
    @test_throws ArgumentError N.burgers_stable_timestep(u,v,-1.,dx,dy)
end

@testset "N07 配布済み熱収支解析" begin
    b=N07Analysis.heat_budget([0.,.25,1.],[2.,4.,5.],[3. 4.;-1. -2.;0. 0.;0. 0.],[.5 .2;0. 0.;-.2 -.4;-.3 -.8])
    @test b.net_input≈[2.,1.]
    @test b.cumulative_input≈[0.,2.,3.]
    @test b.residual≈zeros(3) atol=1e-14
    @test_throws ArgumentError N07Analysis.heat_budget([0.,0.],[1.,1.],zeros(4,1),zeros(4,1))
end
@testset "N07 自作: 開境界2条件の辺別熱収支" begin
    # TODO(自作): 複数ステップで辺別・移流/拡散別の積算、熱量変化、断熱壁を検証する。
    # 入力と独立期待値、許容値の根拠を記す。公式出力は上書きしない。
    @test false
end
@testset "N07 自作: 温度とBurgersの解析解収束" begin
    # TODO(自作): 周期温度の拡散・合成と二成分Burgersを3格子で計算し、誤差と両区間の次数を検証する。
    @test false
end
@testset "N07 保存結果と来歴" begin
    @test check_complete()
end
end
