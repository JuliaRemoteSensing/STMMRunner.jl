module SMUTHI

using CondaPkg
using PythonCall
using ..STMMRunner

const LAYERS = PythonCall.pynew()
const SIMULATION = PythonCall.pynew()
const INITIAL_FIELD = PythonCall.pynew()
const PARTICLES = PythonCall.pynew()
const FAR_FIELD = PythonCall.pynew()
const UTILITY_CUDA = PythonCall.pynew()
const GPU_ENABLED = Ref{Bool}(false)

function enable_gpu()
    if !GPU_ENABLED[]
        CondaPkg.add("cudatoolkit")
        CondaPkg.add("pycuda")
        GPU_ENABLED[] = true
    end
end

function estimate_lmax(x)
    return max(4, Int(ceil(x + 4.05 * ∛x + 2)))
end

function run_smuthi(cfg::STMMConfig)
    if GPU_ENABLED[]
        UTILITY_CUDA.enable_gpu()
    end

    layers = LAYERS.LayerSystem(
        thicknesses=PyList(STMMRunner.thickness.(cfg.layers)),
        refractive_indices=PyList(STMMRunner.refractive_index.(cfg.layers)),
    )

    spheres = PyList([PARTICLES.Sphere(;
        position=[sphere.x, sphere.y, sphere.z],
        refractive_index=sphere.m,
        radius=sphere.r, l_max=estimate_lmax(sphere.r)) for sphere in cfg.spheres])

    plane_wave = INITIAL_FIELD.PlaneWave(
        ; vacuum_wavelength=2π,
        polar_angle=cfg.β * π / 180.0,
        azimuthal_angle=cfg.α * π / 180.0,
        polarization=0
    )

    simulation = SIMULATION.Simulation(;
        layer_system=layers, particle_list=spheres, initial_field=plane_wave)

    simulation.run()

    Csca = pyconvert(Float64, FAR_FIELD.total_scattering_cross_section(;
        simulation=simulation
    ))

    Cext = pyconvert(Float64, FAR_FIELD.extinction_cross_section(;
        simulation=simulation
    ))

    Cabs = Cext - Csca

    return Cext, Csca, Cabs
end

function __init__()
    CondaPkg.add("python"; version="3.9.12")
    CondaPkg.add_pip("smuthi", version="@ git+https://gitlab.com/AmosEgel/smuthi.git")
    PythonCall.pycopy!(LAYERS, pyimport("smuthi.layers"))
    PythonCall.pycopy!(SIMULATION, pyimport("smuthi.simulation"))
    PythonCall.pycopy!(INITIAL_FIELD, pyimport("smuthi.initial_field"))
    PythonCall.pycopy!(PARTICLES, pyimport("smuthi.particles"))
    PythonCall.pycopy!(FAR_FIELD, pyimport("smuthi.postprocessing.far_field"))
    PythonCall.pycopy!(UTILITY_CUDA, pyimport("smuthi.utility.cuda"))
end

end
