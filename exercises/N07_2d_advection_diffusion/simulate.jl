module N07Simulation
using ThermofluidExercise
include("provided_support.jl")
const N=ThermofluidExercise.N07
function simulate_case(id,nx,ny;times=SAVE_TIMES,thermal_step=N.thermal_step!,burgers_step=N.burgers_step!,burgers_timestep=N.burgers_stable_timestep)
    require(nx isa Integer && ny isa Integer && !(nx isa Bool || ny isa Bool) && nx>=3 && ny>=3,"各軸3以上の整数が必要です")
    time=schedule(times); c=case_metadata(id,nx,ny);dx,dy=c["dx"],c["dy"];x=cell_centers(nx,2.);y=cell_centers(ny,1.);isburgers=id=="periodic_burgers"
    a,v=isburgers ? exact_burgers(x,y,0.) : (exact_temperature(x,y,0.,id),nothing)
    b=similar(a);w=isburgers ? similar(v) : nothing; nt=length(time)
    U=Array{Float64}(undef,nx,ny,nt);U[:,:,1]=a
    V=isburgers ? similar(U) : nothing;isburgers && (V[:,:,1]=v)
    adv=zeros(4,nt-1);diff=zeros(4,nt-1);step_time=Float64[];step_dt=Float64[];indices=[0];t=0.
    config=isburgers ? nothing : thermal_config(id)
    dtmax=isburgers ? NaN : N.thermal_stable_timestep(config.cx,config.cy,config.kappa,dx,dy,config.bc)
    if startswith(id,"comparison_")
        ref=thermal_config("periodic_combined");dtmax=N.thermal_stable_timestep(ref.cx,ref.cy,ref.kappa,dx,dy,ref.bc)
    end
    for k in 2:nt
        target=time[k]
        while t<target
            dt=min(isburgers ? burgers_timestep(a,v,.05,dx,dy) : dtmax,target-t)
            require(dt isa Real && isfinite(dt) && dt>0 && t+dt>t,"$id t=$t: 時間が進まないdtです")
            try
                if isburgers
                    burgers_step(b,w,a,v,dt,dx,dy;nu=.05)
                    require(all(isfinite,b) && all(isfinite,w),"非有限速度です")
                    v,w=w,v
                else
                    r=thermal_step(b,a,dt,dx,dy;cx=config.cx,cy=config.cy,kappa=config.kappa,bc=config.bc)
                    require(all(isfinite,b),"非有限温度です")
                    for (j,side) in enumerate(SIDES)
                        adv[j,k-1]+=dt*getproperty(r.boundary_rates.advective,Symbol(side))
                        diff[j,k-1]+=dt*getproperty(r.boundary_rates.diffusive,Symbol(side))
                    end
                end
            catch e
                error("$id $(nx)x$(ny) t=$t dt=$dt: $(sprint(showerror,e))")
            end
            push!(step_time,t);push!(step_dt,dt);t=min(target,t+dt);a,b=b,a
        end
        U[:,:,k]=a;isburgers && (V[:,:,k]=v);push!(indices,length(step_dt))
    end
    merge!(c,Dict("x"=>x,"y"=>y,"time"=>time,"step_time"=>step_time,"step_dt"=>step_dt,"save_step_index"=>indices))
    if isburgers;c["u"]=U;c["v"]=V
    else;c["temperature"]=U;c["advective_heat_integral"]=adv;c["diffusive_heat_integral"]=diff;end
    c
end
function main(;output_dir=DEFAULT_OUTPUT_DIR,thermal_step=N.thermal_step!,burgers_step=N.burgers_step!,burgers_timestep=N.burgers_stable_timestep,publish_options...)
    run_id=bytes2hex(rand(UInt8,16))
    staged(output_dir,("temperature.h5","burgers.h5");publish_options...) do stage
        for family in ("temperature","burgers")
            cases=Dict(id=>Dict(case_id(n...)=>simulate_case(id,n...;thermal_step,burgers_step,burgers_timestep) for n in grids_for(id)) for id in keys(expected_cases(family)))
            write_fields(joinpath(stage,family*".h5"),cases,source_metadata(run_id,family))
        end
        read_pair(stage)
    end
    println("N07の両HDF5を保存しました。analyze.jl、plot.jlを再実行してください。")
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
