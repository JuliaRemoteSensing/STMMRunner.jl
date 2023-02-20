module FaSTMMRunner

export run_fastmm

using DocStringExtensions
using DataFrames: DataFrame
using HDF5: h5open
using UUIDs
using Reexport: @reexport
@reexport using STMMRunner
export run_fastmm

struct FaSTMMOutput
    c_ext::Float64
    c_abs::Float64
    c_scat::Float64
    scattering_matrix::DataFrame
end

"""
$(SIGNATURES)

Use the given configuration to run FaSTMM.

- If `keep = true`, the working directory will not be removed after the run. 
- `fastmm_exe_name` specifies the name or path of your compiled FaSTMM executable.
"""
function run_fastmm(cfg::STMMConfig; keep::Bool = false, fastmm_exe_name::String = "FaSTMM")
    current_dir = pwd()

    id = string(uuid1())
    if cfg.working_directory == ""
        dir = joinpath(pwd(), "FaSTMM_tmp_$id")
    else
        dir = cfg.working_directory
    end
    results = nothing
    proc = nothing

    try
        if !ispath(dir)
            @debug "[Run FaSTMM] Creating workding directory"
            mkdir(dir)
        end

        @debug "[Run FaSTMM] Entering working directory"
        cd(dir)

        @debug "[Run FaSTMM] Writing FaSTMM input file"
        write_geometry(cfg)
        options = make_options(cfg)

        @debug "[Run FaSTMM] Running FaSTMM..."

        # FaSTMM uses OpenMP, so the number of threads needs to be 
        # set via the environment variable.
        ENV["OMP_NUM_THREADS"] = cfg.number_processors
        proc = open(`$fastmm_exe_name $options`, cfg.redirect_stdout; write = true)
        wait(proc)

        @debug "[Run FaSTMM] Collecting FaSTMM output"
        results = collect_output(cfg)
    catch e
        @error "[Run FaSTMM]" exception=(e, catch_backtrace())
        @error "[Run FaSTMM] Runtime error occurred"
    finally
        @debug "[Run FaSTMM] Leaving working directory"
        cd(current_dir)

        if !isnothing(proc)
            @debug "[Run FaSTMM] Killing FaSTMM process"
            kill(proc)
            sleep(0.2)
            if process_running(proc)
                @warn "Failed to kill, please retry manually."
            end
        end

        if !keep && occursin(String(id), dir)
            @assert ispath(dir)
            @debug "[Run FaSTMM] Removing temporary directory"
            rm(dir; force = true, recursive = true)
        end

        return results
    end
end

function write_geometry(cfg::STMMConfig)
    coord = hcat(map(cfg.spheres) do sphere
                     [sphere.x, sphere.y, sphere.z]
                 end...)

    radius = map(s -> s.r, cfg.spheres)
    permittivity = map(s -> s.m^2, cfg.spheres)
    tind = zeros(UInt32, length(radius))
    angles = zeros(length(radius))

    h5open("geometry.h5", "w") do file
        write(file, "/coord", coord)
        write(file, "/radius", radius)
        write(file, "/param_r", real.(permittivity))
        write(file, "/param_i", imag.(permittivity))
        write(file, "/tind", tind)
        write(file, "/angles", angles)
    end
end

function make_options(cfg::STMMConfig)
    options = String[]

    push!(options, "-N_theta")
    push!(options, string(Int(180.0 / cfg.Δθ) + 1))

    push!(options, "-N_phi")
    push!(options, string(Int(360.0 / cfg.Δϕ)))

    push!(options, "-max_iter")
    push!(options, string(cfg.max_iterations))

    push!(options, "-tol")
    push!(options, string(cfg.solution_epsilon))

    if cfg.use_monte_carlo_integration
        push!(options, "-N_ave")
        push!(options, string(cfg.number_incident_directions))
    end

    return options
end

function collect_output(cfg::STMMConfig)
    c_ext = 0.0
    c_abs = 0.0
    c_scat = 0.0

    s = h5open("mueller.h5", "r") do file
        crs = read(file, "/cross_sections")
        c_ext = crs[1]
        c_scat = crs[2]
        c_abs = crs[3]

        mueller = read(file, "/mueller")
        _, c = size(mueller)

        θ = mueller[:, end - 16] * 180.0 / π
        s11 = mueller[:, end - 15]
        s12 = mueller[:, end - 14]
        s13 = mueller[:, end - 13]
        s14 = mueller[:, end - 12]
        s21 = mueller[:, end - 11]
        s22 = mueller[:, end - 10]
        s23 = mueller[:, end - 9]
        s24 = mueller[:, end - 8]
        s31 = mueller[:, end - 7]
        s32 = mueller[:, end - 6]
        s33 = mueller[:, end - 5]
        s34 = mueller[:, end - 4]
        s41 = mueller[:, end - 3]
        s42 = mueller[:, end - 2]
        s43 = mueller[:, end - 1]
        s44 = mueller[:, end]

        if cfg.normalize_scattering_matrix
            s12 ./= s11
            s13 ./= s11
            s14 ./= s11
            s21 ./= s11
            s22 ./= s11
            s23 ./= s11
            s24 ./= s11
            s31 ./= s11
            s32 ./= s11
            s33 ./= s11
            s34 ./= s11
            s41 ./= s11
            s42 ./= s11
            s43 ./= s11
            s44 ./= s11
        end

        if c == 18
            ϕ = mueller[:, 1] * 180.0 / π
            return DataFrame(θ = θ,
                             ϕ = ϕ,
                             s11 = s11,
                             s12 = s12,
                             s13 = s13,
                             s14 = s14,
                             s21 = s21,
                             s22 = s22,
                             s23 = s23,
                             s24 = s24,
                             s31 = s31,
                             s32 = s32,
                             s33 = s33,
                             s34 = s34,
                             s41 = s41,
                             s42 = s42,
                             s43 = s43,
                             s44 = s44)
        else
            @assert c == 17

            return DataFrame(θ = θ,
                             s11 = s11,
                             s12 = s12,
                             s13 = s13,
                             s14 = s14,
                             s21 = s21,
                             s22 = s22,
                             s23 = s23,
                             s24 = s24,
                             s31 = s31,
                             s32 = s32,
                             s33 = s33,
                             s34 = s34,
                             s41 = s41,
                             s42 = s42,
                             s43 = s43,
                             s44 = s44)
        end
    end

    return FaSTMMOutput(c_ext, c_abs, c_scat, s)
end

end
