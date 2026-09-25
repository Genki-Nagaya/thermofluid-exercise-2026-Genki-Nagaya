module N06Simulation
using ThermofluidExercise
include("provided_support.jl")
function simulate_case(nx,ny;times=SAVE_TIMES,cx=1.,cy=.5,safety=.8,t_final=1.,stepper=ThermofluidExercise.N06.advection_step!)
    require(nx isa Integer && ny isa Integer && !(nx isa Bool || ny isa Bool) && nx>=3 && ny>=3,"nx,nyは3以上の整数です")
    time=schedule(times,t_final); dx=2/nx; dy=1/ny
    dt_max=ThermofluidExercise.N06.stable_timestep(cx,cy,dx,dy;safety)
    x=collect((0:nx-1).*dx); y=collect((0:ny-1).*dy)
    a=exact_field(x,y,0.;cx,cy); b=similar(a); U=Array{Float64}(undef,nx,ny,length(time)); U[:,:,1]=a
    dts=Float64[]; counts=Int[]
    for k in 2:length(time)
        steps,dt=ThermofluidExercise.fit_timestep(dt_max,time[k]-time[k-1])
        push!(dts,dt); push!(counts,steps)
        for step in 1:steps
            stepper(b,a,dt,dx,dy;cx,cy)
            all(isfinite,b) || error("非有限値: case=$(case_id(nx,ny)), step=$step, t=$(time[k-1]+step*dt), CFL=$(cx*dt/dx+cy*dt/dy)")
            a,b=b,a
        end
        U[:,:,k]=a
    end
    Dict{String,Any}("nx"=>nx,"ny"=>ny,"dx"=>dx,"dy"=>dy,"dt_max"=>dt_max,"x"=>x,"y"=>y,"time"=>time,"u"=>U,"segment_dt"=>dts,"segment_steps"=>counts)
end
function main(;output_dir=DEFAULT_OUTPUT_DIR,grids=OFFICIAL_GRIDS,times=SAVE_TIMES,stepper=ThermofluidExercise.N06.advection_step!,publish_options...)
    staged(output_dir,(FIELD_NAME,);publish_options...) do stage
        cases=Dict(case_id(nx,ny)=>simulate_case(nx,ny;times,stepper) for (nx,ny) in grids)
        write_fields(joinpath(stage,FIELD_NAME),cases)
    end
    println("N06 fields.h5を保存しました。analyze.jl、plot.jlを順に再実行してください。")
    joinpath(output_dir,FIELD_NAME)
end
abspath(PROGRAM_FILE)==(@__FILE__) && main()
end
