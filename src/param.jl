import Configurations: @option, from_kwargs, from_dict

const Position = NamedTuple{(:x, :y, :z),Tuple{Float64,Float64,Float64}}

@option "sphere" struct SphereConfig
    radius::Float64
    center::Vector{Float64}
    real_ref_index::Union{Float64,Nothing} = nothing
    imag_ref_index::Union{Float64,Nothing} = nothing
    real_chiral_factor::Union{Float64,Nothing} = nothing
    imag_chiral_factor::Union{Float64,Nothing} = nothing
    t_matrix::Union{String,Nothing} = nothing
end

Sphere(; kwargs...) = from_kwargs(SphereConfig; kwargs...)

# This patch is added to support nested parsing.
from_dict(::Type{T}, t::T) where {T} = t

@option "run" struct STMMConfig
    working_directory::String = ""
    number_processors::Int = 1
    sphere_position_file::String = ""
    output_file::String = "mstm.out"
    run_print_file::String = "mstm.run"
    write_sphere_data::Bool = true
    length_scale_factor::Float64 = 1.0
    real_ref_index_scale_factor::Float64 = 1.0
    imag_ref_index_scale_factor::Float64 = 1.0
    real_chiral_factor::Float64 = 0.0
    imag_chiral_factor::Float64 = 0.0
    medium_real_ref_index::Float64 = 1.0
    medium_imag_ref_index::Float64 = 0.0
    medium_real_chiral_factor::Float64 = 0.0
    medium_imag_chiral_factor::Float64 = 0.0
    target_euler_angles_deg::Vector{Float64} = [0.0, 0.0, 0.0]
    mie_epsilon::Float64 = 1e-6
    translation_epsilon::Float64 = 1e-6
    solution_epsilon::Float64 = 1e-8
    max_number_iterations::Int64 = 5000
    store_translation_matrix::Bool = false
    sm_number_processors::Int64 = 1
    near_field_distance::Float64 = 1e8
    iterations_per_correction::Int64 = 20
    θₘᵢₙ::Float64 = 0.0
    θₘₐₓ::Float64 = 180.0
    ϕₘᵢₙ::Float64 = 0.0
    ϕₘₐₓ::Float64 = 0.0
    Δθ::Float64 = 1.0
    frame::String = "target"
    normalize_scattering_matrix::Bool = true
    azimuth_average_scattering_matrix::Bool = false
    gaussian_beam_constant::Float64 = 0.0
    gaussian_beam_focal_point::Array{Float64,1} = [0.0, 0.0, 0.0]
    orientation::String = "fixed"
    incident_azimuth_angle_deg::Float64 = 0.0
    incident_polar_angle_deg::Float64 = 0.0
    calculate_scattering_coefficients::Bool = true
    scattering_coefficient_file::String = ""
    track_iterations::Bool = true
    calculate_near_field::Bool = true
    near_field_plane_coord::Int64 = 1
    near_field_plane_position::Float64 = 0.0
    near_field_plane_vertices::Vector{Float64} = [10.0, 10.0]
    spacial_step_size::Float64 = 0.5
    polarization_angle_deg::Float64 = 0.0
    near_field_output_file::String = "mstm.nf"
    near_field_output_data::Int64 = 2
    plane_wave_epsilon::Float64 = 1e-4
    calculate_t_matrix::Bool = true
    t_matrix_file::String = "mstm.tm"
    t_matrix_convergence_epsilon::Float64 = 1e-7
    spheres::Vector{SphereConfig}
end
