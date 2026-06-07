## Time Axes and Dates

In the example below, we create two time axes. One forward in time, the other backward in time.

```julia
using Drifters

time_unit=[:second,:DateTime][2]
time_period=[:month,:day][2]
ST=Drifters.start_times(3, time_unit=time_unit, direction=:forward, period=time_period)

TA_fd=Drifters.TimeAxis(ST[1],ST[2])
TA_bd=Drifters.TimeAxis(ST[2],ST[1])
```

## Velocity Interpolation

The `dxdt!` etc functions compute the tracked individual velocity. 

```@docs
dxdt!
dxy_dt_replay
dxy_dt_CyclicArray
```

## Setup And Postprocessing 

Convenience functions to initialize a simulation and post-process the output are provided. 

```@docs
postprocess_xy
postprocess_MeshArray
add_lonlat!
```

Basic geography support:

```@docs
gcdist
diff
stproj
stproj_inv
randn_lonlat
interp_to_lonlat
interp_to_xy
nearest_to_xy
```

## Toy Problems

These are used to demonstrate and test the package functionalities:

```@docs
random_flow_field
vortex_flow_field
```

## Read External Files

Trajectory simulated by the [MITgcm](https://mitgcm.readthedocs.io/en/latest/?badge=latest) or observed by the [global drifter program](https://www.aoml.noaa.gov/phod/gdp/index.php) can be read from file using, respectively `.read_flt` (from [.jl](https://gaelforget.github.io/.jl/dev/)) or  `OceanRobots.drifters_hourly_read` (from [OceanRobots.jl](https://gaelforget.github.io/OceanRobots.jl/dev/)).
