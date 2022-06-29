@enum Frame IncidentFrame TargetFrame

@enum Orientation FixedOrientation RandomOrientation

"""
A spherical particle.

$(FIELDS)
"""
Base.@kwdef struct Sphere
    "Radius of the sphere."
    r::Float64

    "x coordinate of the sphere."
    x::Float64

    "y coordinate of the sphere."
    y::Float64

    "z coordinate of the sphere."
    z::Float64

    "Complex refractive index of the sphere. Default is `1.0`."
    m::ComplexF64 = 1.0

    "Complex chiral factor of the sphere. Default is `0.0`."
    β::ComplexF64 = 0.0

    "The T matrix file specifying the scattering properties of the sphere. Default is `\"\"` (empty)."
    t_matrix::String = ""
end

"""
A horizontally infinite medium layer with certain thickness.

$(FIELDS)
"""
Base.@kwdef struct Layer
    "Thickness of the layer. Default is `0.0`."
    d::Float64 = 0.0

    "Complex refractive index of the layer. Default is `1.0`."
    m::ComplexF64 = 1.0

    "Complex chiral factor of the layer. Default is `0.0`."
    β::ComplexF64 = 0.0
end

"""
Configuration for an STMM run.

$(FIELDS)
"""
Base.@kwdef struct STMMConfig
    "The folder for storing intermediate results. Default is `\"\"`, which means a temporary folder will be created inside the current directory and be used as the working directory."
    working_directory::String = ""

    "The number of processes to be initiated. Default is `1`."
    number_processors::Int = 1

    "Name of the output file. Default is `\"mstm.out\"`."
    output_file::String = "mstm.out"

    "Name of the file for storing run information. Default is `\"mstm.run\"`."
    run_print_file::String = "mstm.run"

    "[**MSTM v3**] Rotate the target about the target coordinate frame using a z-y-z specified Euler angle. Default is `(0.0, 0.0, 0.0)`."
    target_euler_angles_deg::NTuple{3,Float64} = (0.0, 0.0, 0.0)

    """
    [**MSTM v3**, **MSTM v4**] Convergence criterion for determining the truncation order for the wave function expansions for each sphere. Default is `1e-6`.

    > Can also be set as a negative integer `-L` to force all sphere expansions to be truncated at `L`-th order.  
    """
    mie_epsilon::Float64 = 1e-6

    "[**MSTM v3**, **MSTM v4**] Convergence criterion for estimating the maximum order of the T matrix for the cluster. Default is `1e-6`."
    translation_epsilon::Float64 = 1e-6

    "[**MSTM v3**, **MSTM v4**] Error criterion for the biconjugate gradient iterative solver. Default is `1e-8`."
    solution_epsilon::Float64 = 1e-8

    "[**MSTM v4**] The maximum truncation degree of the T matrix. This is needed to mostly to keep the code from crashing due to memory allocation limits. Default is `120`."
    max_t_matrix_order::Int = 120

    """
    [**MSTM v3**] Convergence criterion for T matrix solution. Default is `1e-7`.

    > MSTM v4 should also support this option, but it is not yet implemented.
    """
    t_matrix_epsilon::Float64 = 1e-7

    "[**MSTM v3**, **MSTM v4**] The maximum number of iterations attempted to reach a solution. Default is `5000`."
    max_iterations::Int64 = 5000

    """
    [**MSTM v3**, **MSTM v4**] `= false`, the near field translation matrices are calculated during iteration and 
    not stored in memory; `= true`, near field matrices are calculated prior to the solution and stored in memory. 
    Storing the matrices in memory results in some improvement in execution time, yet the amount of memory used
    vs. the resources available is not checked. Use this option at your own peril. Default is `false`.
    """
    store_translation_matrix::Bool = false

    """
    [**MSTM v3**] The number of processors to be used on parallel platform runs to calculate the
    expansion coefficients for the random orientation scattering matrix. The actual number of processors
    used during this step will be the minimum of the total number of processors for the run (the MPI size)
    and sm_number_processors. Default is `1`.
    """
    sm_number_processors::Int64 = 1

    """
    [**MSTM v3**] The threshold for the translation from near-field calculation to far-field calculation.

    - A value larger than the maximum sphere-to-sphere ``kr`` (i.e., a sufficiently large number, say 1e6) causes all spheres to belong to the near field set.
    - A value smaller than the minimum sphere-to-sphere ``kr`` (i.e., zero) causes all spheres to belong to the far field set.
    - A negative value enables the automatic criterion:

    ```math
    kr_{i-j}\\le\\left[\\frac{1}{2}(L_i+L_j)^2\\right]
    ```

    Default is `1e6`, which means all spheres are considered to be within each other's near field.
    """
    near_field_translation_distance::Float64 = 1e6

    """
    [**MSTM v3**] The maximum number of BCGM iterations for a given correction stage. Default is `20`.
    Only works when there are some sphere pairs belonging to the far field set.
    """
    iterations_per_correction::Int64 = 20

    "[**MSTM v3**] The minimum scattering polar angle. Default is `0.0`."
    θ_min::Float64 = 0.0

    "[**MSTM v3**] The maximum scattering polar angle. Default is `0.0`."
    θ_max::Float64 = 180.0

    "[**MSTM v3**] The minimum scattering azimuthal angle. Default is `0.0`."
    ϕ_min::Float64 = 0.0

    "[**MSTM v3**] The maximum scattering azimuthal angle. Default is `0.0`."
    ϕ_max::Float64 = 0.0

    "[**MSTM v3**, **MSTM v4**] The interval of θ for scattering matrix calculation. Default is `1.0`."
    Δθ::Float64 = 1.0

    """
    [**MSTM v3**, **MSTM v4**] This selects the reference coordinate frame upon which the scattering angles
    are based. 

    - `IncidentFrame`, the frame is based on the incident field so that `θ = 0` is the forward scattering direction (i.e., the direction of incident field propagation) and `θ = β, φ = 180` is the target z axis direction. 
    - `TargetFrame`, the frame is based on the target frame so that `θ = 0`` corresponds to the target z axis and `θ = β, φ = α`` is the incident field direction .false. Random orientation calculations force this option.
    """
    frame::Frame = TargetFrame

    "[**MSTM v3**, **MSTM v4**] Whether or not to normalize the scattering matrix. Default is `true`."
    normalize_scattering_matrix::Bool = true

    "[**MSTM v3**, **MSTM v4**] Whether or not to take azimuthal average for fixed orientation calculations. Default is `false`."
    azimuthal_average::Bool = false

    """
    [**MSTM v3**, **MSTM v4**] Dimensionless inverse width CB, at the focal point, of an incident Gaussian
    profile beam. Setting to `0`` selects plane wave incidence. The localized approximation used to
    represent the Gaussian beam is accurate for CB ≤ 0.2. Default is `0.0` (plane wave).
    """
    gaussian_beam_constant::Float64 = 0.0

    """
    [**MSTM v3**, **MSTM v4**] x, y, and z coordinates of the focal point of the Gaussian beam, relative to
    the origin defined by `spheres`. Default is `(0.0, 0.0, 0.0)`.
    """
    gaussian_beam_focal_point::NTuple{3,Float64} = (0.0, 0.0, 0.0)

    """
    [**MSTM v3**, **MSTM v4**] 

    - Use `FixedOrientation` for a fixed orientation calculation.
    - Use `RandomOrientation` for a random orientation calculation.

    Default is `FixedOrientation`.
    """
    orientation::Orientation = FixedOrientation

    "[**MSTM v4**] Use Monte Carlo integration instead of analytical average. Default is `false`."
    use_monte_carlo_integration::Bool = false

    "[**MSTM v4**] Number of directions used for Monte Carlo integration. Works only when `use_monte_carlo_integration = true`. Default is `100`."
    number_incident_directions::Int = 100

    "The azimuthal angle (in degrees) of the incident wave. Default is `0.0`."
    α::Float64 = 0.0

    "The polar angle (in degrees) of the incident wave. Default is `0.0`."
    β::Float64 = 0.0

    "[**MSTM v4**] Whether or not to calculate the scattering matrix. Default is `true`."
    calculate_scattering_matrix::Bool = true

    """
    [**MSTM v3**] Whether or not to calculate the scattering oefficients. Default is `true`.
    When `= false`, the scattering coefficients will be read from the file specified by
    `scattering_coefficient_file`.
    """
    calculate_scattering_coefficients::Bool = true

    """
    [**MSTM v3**] File name for the scattering coefficient file. Default is `\"mstm.sca\"`.

    - When `calculate_scattering_coefficients = false`, the scattering coefficients will be read from this file.
    - When `calculate_scattering_coefficients = true`, the scattering coefficients will be written to this file.
    """
    scattering_coefficient_file::String = "mstm.scat"

    "[**MSTM v4**] Whether or not to use the FFT-based accelerated algorithm. Default is `false`."
    use_fft_translation::Bool = false

    """
    [**MSTM v4**]
    This determines the grid size ``d`` of the cubic lattice used in the algorithm. The
    formula is:

    ```math
    d=\\left(\\frac{V_S}{cN_{s,ext}}\\right)^{1/3}
    ```

    where ``V_S`` is the dimensionless total volume of the external spheres and ``c`` is `cell_volume_fraction`.
    The optimum value is usually where the total number of nodes is close to the total number of external
    spheres NS,ext. Setting `cell_volume_fraction = 0` will automate the setting of d using an estimate of
    the sphere volume fraction. Default is `0.0`. Works only when `use_fft_translation = true`.
    """
    cell_volume_fraction::Float64 = 0.0

    """
    [**MSTM v4**]
    Truncation degree for the node–based expansions of the scattered field. Set
    `node_order = −L`, where `L` is a positive integer, will set the node order to `ceil(d) + L`.
    Default is `-1`. Works only when `use_fft_translation = true`.
    """
    node_order::Int = -1

    """
    [**MSTM v4**] `= 0` will have the scattering matrix printed at a range of polar angles
    (θ), with angular increment given in degrees by `Δθ`. The scattering plane
    coincides with the incident azimuth plane. The limits on θ depend on the number of plane boundaries
    NB. scattering_map_model= 1 will print the scattering matrix in a pair (cosθ > 0, < 0) of 2D
    arrays, with coordinates `kx = sinθcosφ`, `ky = sinθsinφ`. The number of points on each axis is given
    by `scattering_map_dimension`. Default is `0`, and random orientation calculations force `0`.
    """
    scattering_map_model::Int = 0

    """
    [**MSTM v4**] The number of points on each axis. Default is `15`. Works only when `scattering_map_model = 1`.
    """
    scattering_map_dimension::Int = 15

    "[**MSTM v3**] Whether or not to print each iteration's error in the run print file. Default is `true`."
    track_iterations::Bool = true

    "[**MSTM v3**, **MSTM v4**] Whether or not to calculate the near field. Default is `false`."
    calculate_near_field::Bool = false

    "[**MSTM v3**] Coordinate of the near field plane. `1` = y-z, `2` = z-x, `3` = x-y. Default is `1`."
    near_field_plane_coord::Int64 = 1

    """
    [**MSTM v3**] Position of the near field plane. 

    - When `near_field_plane_coord = 1`, this is the x-coordinate.
    - When `near_field_plane_coord = 2`, this is the y-coordinate. 
    - When `near_field_plane_coord = 3`, this is the z-coordinate. 

    Default is `0.0`.
    """
    near_field_plane_position::Float64 = 0.0

    """
    [**MSTM v3**] Rectangle region of the near field plane. 

    - When passing `(ΔX, ΔY)`, the actual region would be `(X_min-ΔX, X_max+ΔX) x (Y_min-ΔY, Y_max+ΔY)`, where `X_min`, `X_max`,`Y_min` and `Y_max` are automatically calculated. 
    - When passing `(X_min, Y_min, X_max, Y_max)`, the region would be used as it is.

    Default is `(0.0, 0.0)`.
    """
    near_field_plane_vertices::Union{NTuple{2,Float64},NTuple{4,Float64}} = (0.0, 0.0)

    "[**MSTM v4**] The minimum x coordinate of the near field calculation region. Default is `0.0`."
    near_field_x_min::Float64 = 0.0

    "[**MSTM v4**] The maximum x coordinate of the near field calculation region. Default is `0.0`."
    near_field_x_max::Float64 = 0.0

    "[**MSTM v4**] The minimum y coordinate of the near field calculation region. Default is `0.0`."
    near_field_y_min::Float64 = 0.0

    "[**MSTM v4**] The maximum y coordinate of the near field calculation region. Default is `0.0`."
    near_field_y_max::Float64 = 0.0

    "[**MSTM v4**] The minimum z coordinate of the near field calculation region. Default is `0.0`."
    near_field_z_min::Float64 = 0.0

    "[**MSTM v4**] The maximum z coordinate of the near field calculation region. Default is `0.0`."
    near_field_z_max::Float64 = 0.0

    "[**MSTM v3**, **MSTM v4**] The spacial step size of the grid points. Default is `0.1`."
    near_field_step_size::Float64 = 0.1

    """
    [**MSTM v4**] Controls where the incident field appears in the calculation results; `= 1` 
    is the standard model, where the external field is = scattered + incident and the field 
    inside the spheres is calculated from the internal field expansions; `≠ 1` has the external
    field due solely to the scattered field, and the field inside the particles is now 
    internal-incident. Default is `1`.
    """
    near_field_calculation_model::Int = 1

    """
    [**MSTM v3**] A specific polarization state of the incident field is needed to calculate the near
    field. The field is taken to be linearly polarized, with a polarization angle of `γ` relative to the k–z
    plane. When `β = 0`, `γ` becomes the azimuth angle of the incident electric field vector relative to the
    cluster coordinate system. Default is `0.0`.
    """
    polarization_angle_deg::Float64 = 0.0

    """
    [**MSTM v4**] Whether or not to enable an accelerated algorithm for calculation of field value.
    Default is `true`.
    """
    store_surface_vector::Bool = true

    """
    [**MSTM v4**] Size of the grid used for re–expansion of fields resulting from plane
    boundary and periodic lattice interactions. Default is `5.0`. Experiment with this 
    if your calculations are taking too long. Only works when `store_surface_vector = true`.
    """
    near_field_expansion_spacing::Float64 = 5.0

    """
    [**MSTM v4**] Expansion truncation order used for the grid re–expansion of the secondary
    field, per the above option. Should scale with near_field_expansion_spacing following a
    Mie–type criteria. The bigger the value, the more accurate the results, but the slower the
    calculations and the greater the memory demands. Default is `10`. Only works when 
    `store_surface_vector = true`.
    """
    near_field_expansion_order::Int = 10

    "[**MSTM v3**, **MSTM v4**] File name for the near field output. Default is `\"mstm.nf\"`."
    near_field_output_file::String = "mstm.nf"

    """
    [**MSTM v3**] The incident field component for the near field calculations -- for either the plane
    wave or Gaussian beam models –- is calculated using a single, regular VSWF expansion centered about
    the beam focal point. The plane wave epsilon is a convergence criterion for this expansion. Default is
    `1e-4`.
    """
    plane_wave_epsilon::Float64 = 1e-4

    """
    [**MSTM v3**] Option to calculate the T matrix or use an exisiting T matrix file. Default is `true`,
    which calculates the T matrix.
    """
    calculate_t_matrix::Bool = true

    """
    [**MSTM v3**] 

    - When `calculate_t_matrix = true`, this is the file name for the T matrix output. 
    - When `calculate_t_matrix = false`, this is the file name for reading the T matrix from.

    Default is `\"mstm.tm\"`.
    """
    t_matrix_file::String = "mstm.tm"

    """
    [**MSTM v3**, **MSTM v4**] Medium layers.

    - For **MSTM v3**, only the refractive index and chiral factor of the first layer will be used. Thickness is ignored.
    - For **MSTM v4**, the first layer locates at `z=0` and extends to `-∞`. The last layer extends to `+∞`. The thickness of the first and the last layer is ignored.
    """
    layers::Vector{Layer} = [Layer()]

    "[**MSTM v3**, **MSTM v4**] Spheres."
    spheres::Vector{Sphere} = Sphere[]

    "[**MSTM v4**] When `periodic = true`, the spheres will be considered to construct a periodic lattice in the x-y plane. Note that periodic lattice is incompatible with random orientation."
    periodic::Bool = false

    "[**MSTM v4**] When `periodic = true`, this option sets the size of the periodic lattice. The cell size will be automatically enlarged if it cannot contain all the spheres."
    cell_size::Tuple{Float64,Float64} = (0.0, 0.0)
end

"""
Helper function for building a configuration on basis of an existing configuration.

$(SIGNATURES)

"""
function STMMConfig(original::STMMConfig; kwargs...)
    dict = Dict(key => getfield(original, key) for key in propertynames(original))
    for (key, value) in kwargs
        dict[key] = value
    end
    return STMMConfig(; dict...)
end

export Orientation, FixedOrientation, RandomOrientation
export Frame, IncidentFrame, TargetFrame
export Sphere, Layer, STMMConfig
