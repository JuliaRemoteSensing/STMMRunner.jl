@testset "Test MSTM v4" begin
    using STMMRunner
    using STMMRunner.MSTM.V4

    sphere1 = Sphere(r = 1.0, x = 4.2, y = 4.3, z = 4.0, m = 1.5)
    sphere2 = Sphere(r = 1.0, x = 0.0, y = 0.3, z = 1.0, m = 1.5)

    spheres = [sphere1, sphere2]
    layers = [Layer(), Layer(d = 10.0, m = 2.0), Layer(m = 1.2)]

    base = STMMConfig(
        number_processors = Sys.CPU_THREADS ÷ 2,
        run_print_file = "",
        working_directory = "",
        scattering_map_model = 0,
        α = 5.0,
        β = 12.0,
        spheres = spheres,
        frame = IncidentFrame,
    )

    @testset "Fixed orientation" begin
        fixed_base = STMMConfig(base;
            orientation = FixedOrientation,
            calculate_near_field = true,
            near_field_x_min = -5.0,
            near_field_x_max = 5.0,
            near_field_y_min = -5.0,
            near_field_y_max = 5.0,
            near_field_z_min = -5.0,
            near_field_z_max = 5.0,
            near_field_step_size = 1.0
        )

        @testset "No layer boundaries, scattering_map_model = 0, frame = incident" begin
            param = STMMConfig(fixed_base)
            @test !isnothing(run_mstm(param))
        end

        @testset "No layer boundaries, scattering_map_model = 0, frame = target" begin
            param = STMMConfig(fixed_base; frame = TargetFrame)
            @test !isnothing(run_mstm(param))
        end

        @testset "No layer boundaries, scattering_map_model = 1, frame = incident" begin
            param = STMMConfig(fixed_base; scattering_map_model = 1, scattering_map_dimension = 20)
            @test !isnothing(run_mstm(param))
        end

        @testset "With layer boundaries, scattering_map_model = 0, frame = incident" begin
            param = STMMConfig(fixed_base; layers = layers)
            @test !isnothing(run_mstm(param))
        end

        @testset "With boundaries, scattering_map_model = 1, frame = incident" begin
            param = STMMConfig(fixed_base; layers = layers, scattering_map_model = 1, scattering_map_dimension = 20)
            @test !isnothing(run_mstm(param))
        end
    end

    @testset "Random orientation" begin
        random_base = STMMConfig(base;
            orientation = RandomOrientation
        )

        @testset "frame = incident" begin
            param = STMMConfig(random_base)
            @test !isnothing(run_mstm(param))
        end

        @testset "frame = target" begin
            param = STMMConfig(random_base; frame = TargetFrame)
            @test !isnothing(run_mstm(param))
        end

        @testset "use Monte Carlo" begin
            param = STMMConfig(random_base; use_monte_carlo_integration = true, number_incident_directions = 10)
            @test !isnothing(run_mstm(param))
        end
    end
end
