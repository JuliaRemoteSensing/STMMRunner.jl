module V3

using DataFrames
using DocStringExtensions
using Printf
using UUIDs
using ...STMMRunner
export run_mstm

abstract type MSTMOutput end

struct NearField
    spheres::DataFrame
    field::DataFrame
end

struct MSTMFixedOutput <: MSTMOutput
    q_ext_unpolarized::Float64
    q_abs_unpolarized::Float64
    q_scat_unpolarized::Float64
    asymmetric_parameter_unpolarized::Float64
    q_ext_parallel::Float64
    q_abs_parallel::Float64
    q_scat_parallel::Float64
    q_ext_perpendicular::Float64
    q_abs_perpendicular::Float64
    q_scat_perpendicular::Float64
    scattering_matrix::DataFrame
    near_field::Union{NearField,Nothing}
end

struct MSTMRandomOutput <: MSTMOutput
    q_ext::Float64
    q_abs::Float64
    q_scat::Float64
    asymmetric_parameter::Float64
    scattering_matrix::DataFrame
    scattering_matrix_expansion_coefficients::DataFrame
    t_matrix::Union{DataFrame,Nothing}
end

"""
Use the given configuration to run MSTM v3.

$(SIGNATURES)

- If `keep = true`, the working directory will not be removed after the run. 
- `mstm_exe_name` specifies the name or path of your compiled MSTM v3 executable.
"""
function run_mstm(cfg::STMMConfig; keep::Bool = false, mstm_exe_name::String = "mstm3")
    current_dir = pwd()

    id = string(uuid1())
    if cfg.working_directory == ""
        dir = joinpath(pwd(), "mstm_tmp_$id")
    else
        dir = cfg.working_directory
    end
    results = nothing
    proc = nothing

    try
        if !ispath(dir)
            @debug "[Run MSTM] Creating workding directory"
            mkdir(dir)
        end

        @debug "[Run MSTM] Entering working directory"
        cd(dir)

        @debug "[Run MSTM] Writing MSTM input file"
        write_input(cfg)

        @debug "[Run MSTM] Running MSTM..."
        proc = open(
            `mpiexec -n $(cfg.number_processors) $(mstm_exe_name)`,
            stdout;
            write = true,
        )
        wait(proc)

        @debug "[Run MSTM] Collecting MSTM output"
        results = collect_output(cfg)
    catch e
        @error "[Run MSTM] $e"
        @error "[Run MSTM] Runtime error occurred"
    finally
        @debug "[Run MSTM] Leaving working directory"
        cd(current_dir)

        if !isnothing(proc)
            @debug "[Run MSTM] Killing MSTM process"
            kill(proc)
            sleep(0.2)
            kill(proc) # We need two SIGINTs to actually kill the process.
            sleep(0.2)
            if process_running(proc)
                @warn "Failed to kill, please retry manually."
            end
        end

        if !keep && occursin(String(id), dir)
            @assert ispath(dir)
            @debug "[Run MSTM] Removing temporary directory"
            rm(dir; force = true, recursive = true)
        end

        return results
    end
end

function format_sphere(sphere::Sphere)::String
    radius_and_center = @sprintf "%.4e %.4e %.4e %.4e " sphere.r sphere.x sphere.y sphere.z

    if !isempty(sphere.t_matrix)
        return radius_and_center * sphere.t_matrix
    else
        return radius_and_center *
               (@sprintf "%.4e %.4e %.4e %.4e" sphere.m.re sphere.m.im sphere.β.re sphere.β.im)
    end
end

# We set `write_sphere_data` to `0` since we do not need them.
# We set `near_field_output_data` to `2` (the full set).
# The manual uses `near_field_translation_distance`, which should be `near_field_distance`.
function write_input(cfg::STMMConfig)
    inp = """
number_spheres
$(length(cfg.spheres))
output_file
$(cfg.output_file)
run_print_file
$(cfg.run_print_file)
write_sphere_data
0
medium_real_ref_index
$(cfg.layers[1].m.re)
medium_imag_ref_index
$(cfg.layers[1].m.im)
medium_real_chiral_factor
$(cfg.layers[1].β.re)
medium_imag_chiral_factor
$(cfg.layers[1].β.im)
target_euler_angles_deg
$(join(cfg.target_euler_angles_deg, ","))
mie_epsilon
$(cfg.mie_epsilon)
translation_epsilon
$(cfg.translation_epsilon)
solution_epsilon
$(cfg.solution_epsilon)
max_number_iterations
$(cfg.max_iterations)
store_translation_matrix
$(Int(cfg.store_translation_matrix))
sm_number_processors
$(cfg.sm_number_processors)
near_field_distance
$(cfg.near_field_translation_distance)
iterations_per_correction
$(cfg.iterations_per_correction)
min_scattering_angle_deg
$(cfg.θ_min)
max_scattering_angle_deg
$(cfg.θ_max)
min_scattering_plane_angle_deg
$(cfg.ϕ_min)
max_scattering_plane_angle_deg
$(cfg.ϕ_max)
delta_scattering_angle_deg
$(cfg.Δθ)
incident_or_target_frame
$(Int(cfg.frame))
normalize_scattering_matrix
$(Int(cfg.normalize_scattering_matrix))
azimuth_average_scattering_matrix
$(Int(cfg.azimuthal_average))
gaussian_beam_constant
$(cfg.gaussian_beam_constant)
gaussian_beam_focal_point
$(join(cfg.gaussian_beam_focal_point, ","))
fixed_or_random_orientation
$(Int(cfg.orientation))
incident_azimuth_angle_deg
$(cfg.α)
incident_polar_angle_deg
$(cfg.β)
calculate_scattering_coefficients
$(Int(cfg.calculate_scattering_coefficients))
scattering_coefficient_file
$(cfg.scattering_coefficient_file)
track_iterations
$(Int(cfg.track_iterations))
calculate_near_field
$(Int(cfg.calculate_near_field))
near_field_plane_coord
$(cfg.near_field_plane_coord)
near_field_plane_position
$(cfg.near_field_plane_position)
near_field_plane_vertices
$(join(cfg.near_field_plane_vertices, ","))
spacial_step_size
$(cfg.near_field_step_size)
polarization_angle_deg
$(cfg.polarization_angle_deg)
near_field_output_file
$(cfg.near_field_output_file)
near_field_output_data
2
plane_wave_epsilon
$(cfg.plane_wave_epsilon)
calculate_t_matrix
$(Int(cfg.calculate_t_matrix))
t_matrix_file
$(cfg.t_matrix_file)
t_matrix_convergence_epsilon
$(cfg.t_matrix_epsilon)
sphere_sizes_and_positions
$(join(map(format_sphere, cfg.spheres),"\n"))
end_of_options
"""

    open("mstm.inp", "w") do io
        write(io, inp)
    end
end

function read_ints(s::AbstractString)::Vector{Int}
    return map(x -> parse(Int, x), filter(!isempty, split(s, r"\s+")))
end

function read_floats(s::AbstractString)::Vector{Float64}
    return map(x -> parse(Float64, x), filter(!isempty, split(s, r"\s+")))
end

function collect_output(cfg::STMMConfig)::MSTMOutput
    @assert isfile(cfg.output_file)

    out = split(read(open(cfg.output_file), String), "\n")
    i = 1

    if cfg.orientation == FixedOrientation
        q = zeros(10)
    else
        q = zeros(4)
    end

    θ = Float64[]
    ϕ = Float64[]
    s11 = Float64[]
    s12 = Float64[]
    s13 = Float64[]
    s14 = Float64[]
    s21 = Float64[]
    s22 = Float64[]
    s23 = Float64[]
    s24 = Float64[]
    s31 = Float64[]
    s32 = Float64[]
    s33 = Float64[]
    s34 = Float64[]
    s41 = Float64[]
    s42 = Float64[]
    s43 = Float64[]
    s44 = Float64[]

    while !occursin("calculation results for run", out[i])
        i += 1
    end

    # Skip header
    i += 2

    if cfg.orientation == FixedOrientation
        # Fixed orientation

        # Read unpolarized efficiencies
        q[1:4] = read_floats(out[i+1])
        i += 2

        # Read parallel efficiencies
        q[5:7] = read_floats(out[i+1])
        i += 2

        # Read perpendicular efficiencies
        q[8:10] = read_floats(out[i+1])
        i += 4

        # Read scattering matrix
        has_phi = !cfg.azimuthal_average
        while !isempty(out[i])
            v = read_floats(out[i])
            i += 1

            push!(θ, v[1])
            push!(ϕ, has_phi ? v[2] : 0.0)

            push!(s44, v[end])
            push!(s43, v[end-1])
            push!(s42, v[end-2])
            push!(s41, v[end-3])
            push!(s34, v[end-4])
            push!(s33, v[end-5])
            push!(s32, v[end-6])
            push!(s31, v[end-7])
            push!(s24, v[end-8])
            push!(s23, v[end-9])
            push!(s22, v[end-10])
            push!(s21, v[end-11])
            push!(s14, v[end-12])
            push!(s13, v[end-13])
            push!(s12, v[end-14])
            push!(s11, v[end-15])
        end
        scattering_matrix = DataFrame(;
            θ,
            ϕ,
            s11,
            s12,
            s13,
            s14,
            s21,
            s22,
            s23,
            s24,
            s31,
            s32,
            s33,
            s34,
            s41,
            s42,
            s43,
            s44,
        )

        if cfg.calculate_near_field && isfile(cfg.near_field_output_file)
            nf = split(read(open(cfg.near_field_output_file), String), "\n")

            # Read near field grid size
            Nx, Ny = read_ints(nf[1])
            Ns = read_ints(nf[2])[1]

            sx = Float64[]
            sy = Float64[]
            sz = fill(cfg.near_field_plane_position, Ns)
            r = Float64[]
            for j = 1:Ns
                v = read_floats(nf[j+2])
                push!(sx, v[1])
                push!(sy, v[2])
                push!(r, v[3])
            end

            if cfg.near_field_plane_coord == 1
                spheres = DataFrame(; x = sz, y = sx, z = sy, r)
            elseif cfg.near_field_plane_coord == 2
                spheres = DataFrame(; x = sy, y = sz, z = sx, r)
            else
                @assert cfg.near_field_plane_coord == 3
                spheres = DataFrame(; x = sx, y = sy, z = sz, r)
            end

            # Read near field data
            E² = Float64[]
            Ex = ComplexF64[]
            Ey = ComplexF64[]
            Ez = ComplexF64[]
            Hx = ComplexF64[]
            Hy = ComplexF64[]
            Hz = ComplexF64[]
            X = Float64[]
            Y = Float64[]
            Z = fill(cfg.near_field_plane_position, Nx * Ny)
            for j = 1:Nx*Ny
                v = read_floats(nf[j+2+Ns])

                push!(X, v[1])
                push!(Y, v[2])

                push!(Ex, complex(v[3], v[4]))
                push!(Ey, complex(v[5], v[6]))
                push!(Ez, complex(v[7], v[8]))

                push!(Hx, complex(v[9], v[10]))
                push!(Hy, complex(v[11], v[12]))
                push!(Hz, complex(v[13], v[14]))
            end

            if cfg.near_field_plane_coord == 1
                x, y, z = Z, X, Y
            elseif cfg.near_field_plane_coord == 2
                x, y, z = Y, Z, X
            else
                @assert cfg.near_field_plane_coord == 3
                x, y, z = X, Y, Z
            end

            field = DataFrame(; x, y, z, Ex, Ey, Ez, Hx, Hy, Hz)
            near_field = NearField(spheres, field)
        else
            near_field = nothing
        end

        return MSTMFixedOutput(q..., scattering_matrix, near_field)
    else
        # Random orientation

        # Read efficiencies
        q[1:4] = read_floats(out[i+1])
        i += 3

        # Read scattering matrix
        Nθ = Int(round((cfg.θ_max - cfg.θ_min) / cfg.Δθ)) + 1
        for j = 1:Nθ
            v = read_floats(out[i+j])

            push!(θ, v[1])
            push!(ϕ, 0.0)

            # Refer to Mishchenko and Yurkin (2017) for the reciprocity relation
            push!(s11, v[2])
            push!(s12, v[3])
            push!(s13, v[4])
            push!(s14, v[5])
            push!(s21, v[3])
            push!(s22, v[6])
            push!(s23, v[7])
            push!(s24, v[8])
            push!(s31, -v[4])
            push!(s32, -v[7])
            push!(s33, v[9])
            push!(s34, v[10])
            push!(s41, v[5])
            push!(s42, v[8])
            push!(s43, -v[10])
            push!(s44, v[11])
        end
        scattering_matrix = DataFrame(;
            θ,
            ϕ,
            s11,
            s12,
            s13,
            s14,
            s21,
            s22,
            s23,
            s24,
            s31,
            s32,
            s33,
            s34,
            s41,
            s42,
            s43,
            s44,
        )
        i += Nθ + 3

        # Read scattering matrix expansion coefficients
        w = Int[]
        a11 = Float64[]
        a12 = Float64[]
        a13 = Float64[]
        a14 = Float64[]
        a22 = Float64[]
        a23 = Float64[]
        a24 = Float64[]
        a32 = Float64[]
        a33 = Float64[]
        a34 = Float64[]
        a44 = Float64[]

        while !isempty(out[i])
            v = read_floats(out[i])
            i += 1

            push!(w, Int(v[1]))
            push!(a11, v[2])
            push!(a22, v[3])
            push!(a33, v[4])
            push!(a23, v[5])
            push!(a32, v[6])
            push!(a44, v[7])
            push!(a12, v[8])
            push!(a34, v[9])
            push!(a13, v[10])
            push!(a24, v[11])
            push!(a14, v[12])
        end
        scattering_matrix_expansion_coefficients =
            DataFrame(; w, a11, a12, a13, a14, a22, a23, a24, a32, a33, a34, a44)

        return MSTMRandomOutput(
            q...,
            scattering_matrix,
            scattering_matrix_expansion_coefficients,
            nothing,
        )
    end
end

end
