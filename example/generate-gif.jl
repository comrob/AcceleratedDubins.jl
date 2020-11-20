#!/usr/bin/env julia

using AcceleratedDubins
using Printf
using PyPlot

for progrun in 1:1
    mkpath("images")
    
    ##############################
    ### Vehicle settings
    r_min = 1.
    r_max = 5.
    params = [5, 15, 2.6, 4]

    ##############################
    ### Position settings

    alpha = 2*pi*rand()
    start = [10*rand(), -5+10*rand(), alpha]
    stop = [10+20*rand(), -5+10*rand(), 0]

    betaarr = 1:0.05:2*pi

    id = 0
    for beta in betaarr
        @printf("%d / %d\n", id, length(betaarr))
        stop[3] = beta

        fig, ax = plt.subplots(1, 2, num="Fastest vs. shortest")
        fig.set_size_inches(6.4*2, 4.8)
        ax[1].axis("equal")
        ax[1].set_xlabel("x")
        ax[1].set_ylabel("y")
        ax[2].set_ylim(-7, 7)


        ax[2].axis("equal")
        ax[2].set_ylabel("Speed [m/s]")
        ax[2].set_xlabel("Time [s]")
        ax[2].set_ylim(params[1], params[2])
        ax[2].set_xlim(-0.5, 5)

        # simple path computation
        # compute fastest path
        # tries all radii combinations
        path_fastest, _ = fastest_path(start, stop, radii_samples_exp(r_min, r_max, 5), params)
        # tries optimized path
        path_optimized, _ = optimized_path(start, stop, [r_min, r_max], params)
        # compute shortest path
        path_shortest, _ = retrieve_path(start, stop, [r_min, r_min], params)
        # sample the path for plotting
        confx, confy = sample_path(path_fastest)
        ax[1].plot(confx, confy, color="green", label="fastest")
        confx, confy = sample_path(path_optimized)
        ax[1].plot(confx, confy, color="blue", label="optimized")
        confx, confy = sample_path(path_shortest)
        ax[1].plot(confx, confy, color="orange", label="shortest")
        

        times, speeds = speed_profile(path_fastest, params)
        ax[2].plot(times, speeds, color="green", label="fastest")
        times, speeds = speed_profile(path_optimized, params)
        ax[2].plot(times, speeds, color="blue", label="optimized")
        times, speeds = speed_profile(path_shortest, params)
        ax[2].plot(times, speeds, color="orange", label="shortest")

        ax[2].legend(bbox_to_anchor=(1.3, 0.999))

        time_fastest = path_time(path_fastest, params)
        time_fastest_no_opt = (sum(path_fastest.lengths) / speed_by_radius(r_min, params[1], params[2]))

        time_optimized = path_time(path_optimized, params)
        time_optimized_no_opt = (sum(path_optimized.lengths) / speed_by_radius(r_min, params[1], params[2]))

        time_shortest = path_time(path_shortest, params)
        time_shortest_no_opt = (sum(path_shortest.lengths) / speed_by_radius(r_min, params[1], params[2]))

        fig.savefig(@sprintf("images/img%03d.jpg", id), bbox_inches=nothing, pad_inches=-0.1)
        id += 1
        # println("[GREEN] Fastest path:    vp = ", round(time_fastest, digits=3), "    no_opt = ", round(time_fastest_no_opt, digits=3))
        # println("[BLUE] Optimized path:   vp = ", round(time_optimized, digits=3), "    no_opt = ", round(time_optimized_no_opt, digits=3))
        # println("[ORANGE] Shortest path:  vp = ", round(time_shortest, digits=3), "    no_opt = ", round(time_shortest_no_opt, digits=3))
        fig.clear()
    end

    gifconvert = `convert -delay 12 -loop 0 "images/*.jpg" "images/myimage.gif"`
    @show gifconvert
    run(gifconvert);
    @printf("Converted.")

end
