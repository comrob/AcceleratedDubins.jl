using Test
using AcceleratedDubins
using Dubins

for run in 1:1
    tol = 1e-6    # for isapprox comparison
    params = [5, 15, 2.6, 4]

    @testset "original vs. accelerated" begin
        # compare lengths of the original path (pkg Dubins) and 
        # accelerated with minimum turning radius (pkg AcceleratedDubins)
        start = [0, 0, 0.5]
        stop = [3, 0, pi]
        r_min = 1.

        # same positions, but shifted and rotated (check is rotation works correctly)
        start_r = [5, 5, 0.5+pi/2] 
        stop_r = [5, 8, pi+pi/2]

        _, path_d = dubins_shortest_path(start, stop, r_min)
        path_a, _ = retrieve_path(start, stop, [r_min, r_min], params) # shortest path
        
        _, path_d_r = dubins_shortest_path(start_r, stop_r, r_min)
        path_a_r, _ = retrieve_path(start_r, stop_r, [r_min, r_min], params) # shortest path
        

        @test isapprox(dubins_path_length(path_d), path_len(path_a), atol=tol)
        @test isapprox(dubins_path_length(path_d_r), path_len(path_a), atol=tol)

    end

    @testset "straight path" begin
        start = [0, 0, 0.]
        stop = [6, 0, 0.]
        
        _, path_d = dubins_shortest_path(start, stop, 1.)
        len_d = dubins_path_length(path_d)

        path_s, _ = retrieve_path(start, stop, [1., 1.], params) # shortest path
        path_f, _ = fastest_path(start, stop, [1., 2., 3.], params)
        path_o, _ = optimized_path(start, stop, [1., 3.], params) 

        len_s = path_len(path_s)
        len_f = path_len(path_f)
        len_o = path_len(path_o)

        @test isapprox(len_d, path_len(path_s), atol=tol)
        @test isapprox(len_s, len_f, atol=tol)
        @test isapprox(len_s, len_o, atol=tol)
    end

    @testset "proposed paths better or equal" begin
        for testrun in 1:10
            start = [10*rand(), 10*rand(), 2*pi*rand()]
            stop = [10*rand()+10, 10*rand()+5, 2*pi*rand()]
            
            _, path_d = dubins_shortest_path(start, stop, 1.)
            len_d = dubins_path_length(path_d)

            path_s, _ = retrieve_path(start, stop, [1., 1.], params) # shortest path
            path_f, _ = fastest_path(start, stop, [1., 2., 3.], params)
            path_o, _ = optimized_path(start, stop, [1., 3.], params) 

            time_s = path_time(path_s, params)
            time_f = path_time(path_f, params)
            time_o = path_time(path_o, params)

            @test time_s >= time_f
            @test time_s >= time_o
        end
    end
end
