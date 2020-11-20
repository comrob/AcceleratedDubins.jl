using AcceleratedDubins

## Parameters of the vehicle
params = [5., 15., 2.6, 4] ## v_min v_max a_max -a_min
r_min = 1. # minimum turning radius
r_max = 5. # maximum turning radius

## Generate random configurations
start = [10*rand(), -10+5*rand(), 2*pi*rand()]
stop = [10+20*rand(), -10+5*rand(), 2*pi*rand()]

## try radii combinations and get resulting path time
    path_fastest, _ = fastest_path(start, stop, radii_samples_exp(r_min, r_max, 3), params)
    time_fastest = path_time(path_fastest, params)
## try to optimize path using Optim
    path_optimized, _ = optimized_path(start, stop, [r_min, r_max], params)
    time_optimized = path_time(path_optimized, params)
## compute shortest path (original Dubins with r_min)
    path_shortest, _ = retrieve_path(start, stop, [r_min, r_min], params)
    time_shortest = path_time(path_shortest, params)

@show time_fastest, time_optimized, time_shortest

## sample the path for plotting 
times, speeds = speed_profile(path_fastest, params)
