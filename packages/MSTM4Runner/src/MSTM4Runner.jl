module MSTM4Runner

using MSTM_jll: mpiexec, mstm
using DataFrames
using DocStringExtensions
using Printf
using Scanf
using UUIDs
using Reexport: @reexport
@reexport using STMMRunner
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
    q_ext_parallel::Float64
    q_abs_parallel::Float64
    q_scat_parallel::Float64
    q_ext_perpendicular::Float64
    q_abs_perpendicular::Float64
    q_scat_perpendicular::Float64
    q_scat_down_unpolarized::Float64
    q_scat_up_unpolarized::Float64
    q_scat_down_parallel::Float64
    q_scat_up_parallel::Float64
    q_scat_down_perpendicular::Float64
    q_scat_up_perpendicular::Float64
    q_ext_down_unpolarized::Float64
    q_ext_up_unpolarized::Float64
    q_ext_down_parallel::Float64
    q_ext_up_parallel::Float64
    q_ext_down_perpendicular::Float64
    q_ext_up_perpendicular::Float64
    q_scat_waveguide_unpolarized::Float64
    q_scat_waveguide_parallel::Float64
    q_scat_waveguide_perpendicular::Float64
    scattering_matrix::DataFrame
    near_field::Union{NearField, Nothing}
end

struct MSTMLatticeOutput <: MSTMOutput
    ref_unpolarized::Float64
    abs_unpolarized::Float64
    scat_unpolarized::Float64
    ref_parallel::Float64
    abs_parallel::Float64
    scat_parallel::Float64
    ref_perpendicular::Float64
    abs_perpendicular::Float64
    scat_perpendicular::Float64
    scattering_matrix::DataFrame
end

struct MSTMRandomOutput <: MSTMOutput
    q_ext::Float64
    q_abs::Float64
    q_scat::Float64
    scattering_matrix::DataFrame
end

struct MSTMMonteCarloOutput <: MSTMOutput
    q_ext::Float64
    q_abs::Float64
    q_scat::Float64
    q_scat_down_unpolarized::Float64
    q_scat_up_unpolarized::Float64
    q_scat_down_parallel::Float64
    q_scat_up_parallel::Float64
    q_scat_down_perpendicular::Float64
    q_scat_up_perpendicular::Float64
    scattering_matrix_azimuthal_average::DataFrame
    scattering_matrix_diffuse::DataFrame
    scattering_matrix_expansion_coefficients::DataFrame
end

"""
$(SIGNATURES)

Use the given configuration to run MSTM v4.

- If `keep = true`, the working directory will not be removed after the run. 
"""
function run_mstm(cfg::STMMConfig; keep::Bool = false,
                  mstm_command::Union{Cmd, Nothing} = nothing)
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

        # We cannot simply do `$(mpiexec()) ... $(mstm())` here due to the limit of Cmd interpolation.
        # See https://github.com/JuliaLang/julia/issues/39282 for more details.
        @debug "[Run MSTM] Running MSTM..."
        proc = open(isnothing(mstm_command) ?
                    setenv(`$(mpiexec().exec[1]) -n $(cfg.number_processors) $(mstm().exec[1])`,
                           mstm().env) : mstm_command,
                    cfg.redirect_stdout;
                    write = true)
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
    radius_and_center = @sprintf "%.4e %.4e %.4e %.4e " sphere.x sphere.y sphere.z sphere.r

    if !isempty(sphere.t_matrix)
        return radius_and_center * sphere.t_matrix
    else
        mₗ = sphere.m / (1.0 - sphere.β * sphere.m)
        mᵣ = sphere.m / (1.0 + sphere.β * sphere.m)

        return radius_and_center *
               (@sprintf "(%.4e,%.4e) (%.4e,%.4e)" mₗ.re mₗ.im mᵣ.re mᵣ.im)
    end
end

function format_spheres(cfg::STMMConfig)
    return """
number_spheres
$(length(cfg.spheres))
sphere_data
$(join(map(format_sphere, cfg.spheres),"\n"))
end_of_sphere_data"""
end

function format_fft_options(cfg::STMMConfig)
    if !cfg.use_fft_translation
        return "fft_translation_option\nfalse"
    else
        return """
        fft_translation_option
        true
        min_fft_nsphere
        $(length(cfg.spheres))
        cell_volume_fraction
        $(cfg.cell_volume_fraction)
        node_order
        $(cfg.node_order)"""
    end
end

# FIXME Optically active layers are not supported in MSTM v4 yet.
# See https://github.com/dmckwski/MSTM/issues/6
function format_layers(cfg::STMMConfig)::String
    N = length(cfg.layers)

    return """
number_plane_boundaries
$(N - 1)
layer_thickness
$(N >= 3 ? join(map(x -> x.d, cfg.layers[2:N-1]), ",") : "")
layer_ref_index
$(join(map(x -> (@sprintf "(%.4e,%.4e)" x.m.re x.m.im), cfg.layers), ","))"""
end

function format_orientation(cfg::STMMConfig)::String
    if cfg.orientation == RandomOrientation
        if cfg.use_monte_carlo_integration
            return "incidence_average\ntrue\nnumber_incident_directions\n$(cfg.number_incident_directions)"
        else
            return "random_orientation\ntrue"
        end
    else
        return "incident_alpha_deg\n$(cfg.α)\nincident_beta_deg\n$(cfg.β)\nazimuthal_average\n$(cfg.azimuthal_average)\nsingle_origin_expansion\n$(cfg.single_origin_expansion)"
    end
end

function format_near_field(cfg::STMMConfig)::String
    if cfg.calculate_near_field && cfg.orientation == FixedOrientation
        return """
    calculate_near_field
    true
    near_field_output_file
    $(cfg.near_field_output_file)
    near_field_minimum_border
    $(cfg.near_field_x_min),$(cfg.near_field_y_min),$(cfg.near_field_z_min)
    near_field_maximum_border
    $(cfg.near_field_x_max),$(cfg.near_field_y_max),$(cfg.near_field_z_max)
    near_field_step_size
    $(cfg.near_field_step_size)
    near_field_calculation_model
    $(cfg.near_field_calculation_model)
    store_surface_vector
    $(cfg.store_surface_vector)""" * (cfg.store_surface_vector ? """
    \nnear_field_expansion_spacing
    $(cfg.near_field_expansion_spacing)
    near_field_expansion_order
    $(cfg.near_field_expansion_order)""" :
                "")
    else
        return "calculate_near_field\nfalse"
    end
end

function format_lattice(cfg::STMMConfig)
    if !cfg.periodic
        return "periodic_lattice\nfalse"
    else
        if cfg.orientation == RandomOrientation
            throw(ErrorException("Lattice mode is incompatible with random operation!"))
        end

        x_min = minimum(s.x - s.r for s in cfg.spheres)
        x_max = maximum(s.x + s.r for s in cfg.spheres)
        y_min = minimum(s.y - s.r for s in cfg.spheres)
        y_max = maximum(s.y + s.r for s in cfg.spheres)
        dx = x_max - x_min
        dy = y_max - y_min
        cx = max(dx, cfg.cell_size[1])
        cy = max(dy, cfg.cell_size[2])
        return "periodic_lattice\ntrue\ncell_width\n$cx,$cy"
    end
end

# We set `print_sphere_data` to `false` since we do not need them.
function write_input(cfg::STMMConfig)
    inp = """
    $(format_spheres(cfg))
    $(format_layers(cfg))
    output_file
    $(cfg.output_file)
    run_file
    $(cfg.run_print_file)
    print_sphere_data
    false
    mie_epsilon
    $(cfg.mie_epsilon)
    translation_epsilon
    $(cfg.translation_epsilon)
    solution_epsilon
    $(cfg.solution_epsilon)
    t_matrix_convergence_epsilon
    $(cfg.t_matrix_epsilon)
    max_iterations
    $(cfg.max_iterations)
    max_t_matrix_order
    $(cfg.max_t_matrix_order)
    store_translation_matrix
    $(cfg.store_translation_matrix)
    $(format_fft_options(cfg))
    $(format_orientation(cfg))
    $(format_near_field(cfg))
    incident_frame
    $(cfg.frame == IncidentFrame)
    normalize_s11
    $(cfg.normalize_scattering_matrix)
    gaussian_beam_constant
    $(cfg.gaussian_beam_constant)
    gaussian_beam_focal_point
    $(@sprintf("(%.4e,%.4e,%.4e)", cfg.gaussian_beam_focal_point...))
    calculate_scattering_matrix
    $(cfg.calculate_scattering_matrix)
    scattering_map_model
    $(cfg.scattering_map_model)
    scattering_map_increment
    $(cfg.Δθ)
    scattering_map_dimension
    $(cfg.scattering_map_dimension)
    $(format_lattice(cfg))
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

    out = readlines(cfg.output_file)
    i = 1

    if cfg.orientation == FixedOrientation
        q = cfg.periodic ? zeros(9) : zeros(24)
    else
        q = cfg.use_monte_carlo_integration ? zeros(9) : zeros(3)
    end

    θ = Float64[]
    kx = Float64[]
    ky = Float64[]
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
    type = String[]

    while !occursin("calculation results for run", out[i])
        i += 1
    end

    # Skip header
    i += cfg.orientation == FixedOrientation ? 4 : 6

    if cfg.orientation == FixedOrientation
        if cfg.periodic
            # Lattice mode
            q[1:9] = collect(@scanf out[i + 1] "%e%e%e%e%e%e%e%e%e" zeros(9)...)[2:end]

            i += 6
        elseif length(cfg.layers) == 1
            # One layer

            # Read total efficiencies
            q[1:9] = read_floats(out[i + 1])
            i += 2

            # Read hemispherical scattering efficiencies
            q[10:15] = read_floats(out[i + 1])
            i += 3
        else
            # Multiple layers

            # Read total efficiencies
            q[1:9] = read_floats(out[i + 1])
            i += 2

            # Read down/up extinction efficiencies
            q[16:21] = read_floats(out[i + 1])
            i += 2

            # Read hemispherical scattering efficiencies
            q[10:15] = read_floats(out[i + 1])
            i += 2

            # Read waveguide scattering efficiencies
            if occursin("waveguide", out[i])
                q[22:24] = read_floats(out[i + 1])
                i += 3
            else
                i += 1
            end
        end

        current_type = "total"
        if occursin("reflection", out[i])
            current_type = "reflection"
            i += 4
        elseif occursin("kx", out[i])
            current_type = "backward"
            i += 1
        else
            i += 3
        end

        # Read scattering matrix
        while i <= length(out)
            if occursin("transmission", out[i])
                current_type = "transmission"
                i += 2
                continue
            elseif occursin("forward", out[i])
                current_type = "forward"
                i += 2
                if cfg.periodic
                    i += 2
                end
                continue
            end

            v = read_floats(out[i])
            i += 1

            push!(type, current_type)

            # TODO: need to check whether lattice output can support theta instead kx-ky convention
            if cfg.scattering_map_model == 0 && !cfg.periodic
                push!(θ, v[1])
            else
                push!(kx, v[1])
                push!(ky, v[2])
            end

            push!(s44, v[end])
            push!(s43, v[end - 1])
            push!(s42, v[end - 2])
            push!(s41, v[end - 3])
            push!(s34, v[end - 4])
            push!(s33, v[end - 5])
            push!(s32, v[end - 6])
            push!(s31, v[end - 7])
            push!(s24, v[end - 8])
            push!(s23, v[end - 9])
            push!(s22, v[end - 10])
            push!(s21, v[end - 11])
            push!(s14, v[end - 12])
            push!(s13, v[end - 13])
            push!(s12, v[end - 14])
            push!(s11, v[end - 15])
        end

        if cfg.scattering_map_model == 0 && !cfg.periodic
            scattering_matrix = DataFrame(;
                                          θ,
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
                                          type)
        else
            scattering_matrix = DataFrame(;
                                          kx,
                                          ky,
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
                                          type)
        end

        if cfg.calculate_near_field && isfile(cfg.near_field_output_file)
            nf = split(read(open(cfg.near_field_output_file), String), "\n")

            Ns = read_ints(nf[3])[1]
            sx = Float64[]
            sy = Float64[]
            sz = Float64[]
            r = Float64[]
            for j in 1:Ns
                v = read_floats(nf[j + 3])
                push!(sx, v[1])
                push!(sy, v[2])
                push!(sz, v[3])
                push!(r, v[4])
            end
            spheres = DataFrame(; x = sx, y = sy, z = sz, r)

            Nb = read_ints(nf[Ns + 4])[1]
            Nx, Ny, Nz = read_ints(nf[Ns + Nb + 7])

            # Read near field data
            x = Float64[]
            y = Float64[]
            z = Float64[]
            Ex₌ = ComplexF64[]
            Ex⊥ = ComplexF64[]
            Ey₌ = ComplexF64[]
            Ey⊥ = ComplexF64[]
            Ez₌ = ComplexF64[]
            Ez⊥ = ComplexF64[]
            Hx₌ = ComplexF64[]
            Hx⊥ = ComplexF64[]
            Hy₌ = ComplexF64[]
            Hy⊥ = ComplexF64[]
            Hz₌ = ComplexF64[]
            Hz⊥ = ComplexF64[]
            for j in 1:(Nx * Ny * Nz)
                v = read_floats(nf[j + 7 + Ns + Nb])
                push!(x, v[1])
                push!(y, v[2])
                push!(z, v[3])
                push!(Ex₌, complex(v[4], v[5]))
                push!(Ey₌, complex(v[6], v[7]))
                push!(Ez₌, complex(v[8], v[9]))
                push!(Hx₌, complex(v[10], v[11]))
                push!(Hy₌, complex(v[12], v[13]))
                push!(Hz₌, complex(v[14], v[15]))
                push!(Ex⊥, complex(v[16], v[17]))
                push!(Ey⊥, complex(v[18], v[19]))
                push!(Ez⊥, complex(v[20], v[21]))
                push!(Hx⊥, complex(v[22], v[23]))
                push!(Hy⊥, complex(v[24], v[25]))
                push!(Hz⊥, complex(v[26], v[27]))
            end

            field = DataFrame(;
                              x,
                              y,
                              z,
                              Ex₌,
                              Ex⊥,
                              Ey₌,
                              Ey⊥,
                              Ez₌,
                              Ez⊥,
                              Hx₌,
                              Hx⊥,
                              Hy₌,
                              Hy⊥,
                              Hz₌,
                              Hz⊥)

            near_field = NearField(spheres, field)
        else
            near_field = nothing
        end

        return cfg.periodic ? MSTMLatticeOutput(q..., scattering_matrix) :
               MSTMFixedOutput(q..., scattering_matrix, near_field)
    elseif !cfg.use_monte_carlo_integration
        # Random orientation, analytical average

        # Read efficiencies
        q[1:3] = read_floats(out[i + 1])
        i += 4

        # Read scattering matrix
        Nθ, sm = read_ints(out[i])
        i += 1

        if sm == 6
            for j in 1:Nθ
                v = read_floats(out[i + j])

                push!(θ, v[1])
                push!(s11, v[2])
                push!(s12, v[3])
                push!(s22, v[4])
                push!(s33, v[5])
                push!(s34, v[6])
                push!(s44, v[7])
            end
            scattering_matrix = DataFrame(; θ, s11, s12, s22, s33, s34, s44)
        else
            for j in 1:Nθ
                v = read_floats(out[i + j])

                push!(θ, v[1])
                push!(s11, v[2])
                push!(s12, v[3])
                push!(s13, v[4])
                push!(s14, v[5])
                push!(s21, v[6])
                push!(s22, v[7])
                push!(s23, v[8])
                push!(s24, v[9])
                push!(s31, v[10])
                push!(s32, v[11])
                push!(s33, v[12])
                push!(s34, v[13])
                push!(s41, v[14])
                push!(s42, v[15])
                push!(s43, v[16])
                push!(s44, v[17])
            end
            scattering_matrix = DataFrame(; θ, s11, s12, s13, s14, s21, s22, s23, s24, s31,
                                          s32, s33, s34, s41, s42, s43, s44)
        end

        return MSTMRandomOutput(q..., scattering_matrix)
    else
        # Random orientation, Monte Carlo

        # Read efficiencies
        q[1:3] = read_floats(out[i + 1])
        i += 2

        q[4:9] = read_floats(out[i + 1])
        i += 4

        # Read scattering matrix
        Nθ, _ = read_ints(out[i])
        i += 1

        for j in 1:Nθ
            v = read_floats(out[i + j])
            push!(θ, v[1])
            push!(s11, v[2])
            push!(s12, v[3])
            push!(s22, v[4])
            push!(s33, v[5])
            push!(s34, v[6])
            push!(s44, v[7])
        end

        scattering_matrix_azimuthal_average = DataFrame(; θ, s11, s12, s22, s33, s34, s44)

        for v in [θ, s11, s12, s22, s33, s34, s44]
            resize!(v, 0)
        end

        i += Nθ + 1
        @assert occursin("diffuse", out[i])
        i += 1

        for j in 1:Nθ
            v = read_floats(out[i + j])
            push!(θ, v[1])
            push!(s11, v[2])
            push!(s12, v[3])
            push!(s22, v[4])
            push!(s33, v[5])
            push!(s34, v[6])
            push!(s44, v[7])
        end

        scattering_matrix_diffuse = DataFrame(; θ, s11, s12, s22, s33, s34, s44)

        i += Nθ + 1
        @assert occursin("azimuthal averaged scattering matrix expansion coefficients",
                         out[i])
        i += 2

        n = Int[]
        a11 = Float64[]
        a44 = Float64[]
        a12 = Float64[]
        a34 = Float64[]
        a22p = Float64[]
        a22m = Float64[]

        while i <= length(out)
            v = read_floats(out[i])

            push!(n, v[1])
            push!(a11, v[2])
            push!(a44, v[3])
            push!(a12, v[4])
            push!(a34, v[5])
            push!(a22p, v[6])
            push!(a22m, v[7])

            i += 1
        end

        scattering_matrix_expansion_coefficients = DataFrame(; n, a11, a44, a12, a34, a22p,
                                                             a22m)

        return MSTMMonteCarloOutput(q...,
                                    scattering_matrix_azimuthal_average,
                                    scattering_matrix_diffuse,
                                    scattering_matrix_expansion_coefficients)
    end
end

end
