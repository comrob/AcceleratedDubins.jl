#!/usr/bin/env julia

include("../src/AcceleratedDubins.jl")
using .AcceleratedDubins

using PyPlot
using PyCall

# Save figures into a directory "images"
save_figures = false
if save_figures
    mkpath("images")
end

for run in 1:1
    println("COMROB: Dubins path with multiple radii demo")
    println("Library examples")


    ##############################
    ### Vehicle settings
    minr = 1.
    maxr = 5.

    vr = 5.
    vs = 15.

    speed_params = [vr, vs, 2.6, 4]


    ##############################
    ### Position settings

    ### tested settings with visible differences
    # alpha = 0.8*pi
    # beta = 0.1*pi
    # start = [6, -7, alpha]
    # stop = [19, -5, beta]

    ### random settings
    alpha = 2*pi*rand()
    beta = 2*pi*rand()
    start = [10*rand(), -10+5*rand(), alpha]
    stop = [10+20*rand(), -10+5*rand(), beta]


    while true
        fig, ax = plt.subplots(1, 2, num="Fastest vs. shortest")
        # fig.set_size_inches(6.4*2, 4.8)
        # plt.clf()
        ax[1].axis("equal")
        ax[1].set_xlabel("x")
        ax[1].set_ylabel("y")

        ax[2].axis("equal")
        ax[2].set_ylabel("Speed [m/s]")
        ax[2].set_xlabel("Time [s]")
        # ax[2].set_ylim(vr, vs)
        # simple path computation
        # compute fastest path
        # tries all radii combinations
        _, path_fastest = fastest_path(start, stop, Array(minr:2.0:maxr), speed_params)
        # tries optimized path
        _, path_optimized = optimized_path(start, stop, minr, maxr, speed_params)
        # compute shortest path
        _, path_shortest = retrieve_path(start, stop, minr, minr, speed_params)
        # sample the path for plotting
        confx, confy = sample_path(path_fastest)
        ax[1].plot(confx, confy, color="green")
        confx, confy = sample_path(path_optimized)
        ax[1].plot(confx, confy, color="blue")
        confx, confy = sample_path(path_shortest)
        ax[1].plot(confx, confy, color="orange")

        times, speeds = speed_profile(path_fastest, speed_params)
        ax[2].plot(times, speeds, color="green")
        times, speeds = speed_profile(path_optimized, speed_params)
        ax[2].plot(times, speeds, color="blue")
        times, speeds = speed_profile(path_shortest, speed_params)
        ax[2].plot(times, speeds, color="orange")

        time_fastest = path_time(path_fastest, speed_params)
        time_fastest_no_opt = (sum(path_fastest.lengths) / speed_by_radius(minr, speed_params[1], speed_params[2]))

        time_optimized = path_time(path_optimized, speed_params)
        time_optimized_no_opt = (sum(path_optimized.lengths) / speed_by_radius(minr, speed_params[1], speed_params[2]))

        time_shortest = path_time(path_shortest, speed_params)
        time_shortest_no_opt = (sum(path_shortest.lengths) / speed_by_radius(minr, speed_params[1], speed_params[2]))

        println("[GREEN] Fastest path:    vp = ", round(time_fastest, digits=3), "    no_opt = ", round(time_fastest_no_opt, digits=3))
        println("[BLUE] Optimized path:   vp = ", round(time_optimized, digits=3), "    no_opt = ", round(time_optimized_no_opt, digits=3))
        println("[ORANGE] Shortest path:  vp = ", round(time_shortest, digits=3), "    no_opt = ", round(time_shortest_no_opt, digits=3))

        plt.pause(0.1)
        println("Press enter for next random demo, q and enter to quit")
        str = readline()
        close()
        if occursin("q", str) || occursin("Q", str)
            break
        end

        alpha = 2*pi*rand()
        beta = 2*pi*rand()
        start = [10*rand(), -10+5*rand(), alpha]
        stop = [10+20*rand(), -10+5*rand(), beta]
    end


end
