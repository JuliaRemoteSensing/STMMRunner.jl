### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ 93737264-8d97-11ec-2d60-336b353c2e54
begin
    import Pkg
    Pkg.activate("../..")

    import GMT

    using DelimitedFiles
    using STMMRunner
    using STMMRunner.MSTM.V3
end

# ╔═╡ 9e0a299f-9c69-4cb1-8de7-62116eb4e404
begin
    spheres = [
        Sphere(r = 0.7 + 0.5 * rand(), x = i, y = j, z = k, m = 1.5 + 1.0im * rand())
        for i = -6.0:3.0:6.0 for j = -6.0:3.0:6.0 for k = -6.0:3.0:6.0
    ]

    param = STMMConfig(
        number_processors = Sys.CPU_THREADS ÷ 2,
        sm_number_processors = Sys.CPU_THREADS ÷ 2,
        run_print_file = "",
        working_directory = "",
        orientation = FixedOrientation,
        α = 70.0,
        calculate_near_field = true,
        near_field_plane_position = 0.0,
        near_field_plane_vertices = (5.0, 5.0),
        near_field_plane_coord = 1,
        near_field_step_size = 0.1,
        spheres = spheres,
    )
end

# ╔═╡ 4f69842a-3658-4345-a6ce-8dcd0c20c564
results = run_mstm(param; keep = false)

# ╔═╡ a826a137-7c2c-4c70-9217-4797ca95b593
let
    nf = results.near_field.field
    y = nf.y |> Set |> collect |> sort
    z = nf.z |> Set |> collect |> sort
    sp = results.near_field.spheres

    grd = GMT.xyz2grd(
        [nf.y nf.z real.(nf.Ex)],
        region = (first(y), last(y), first(z), last(z)),
        inc = "$(y[2]-y[1])/$(z[2]-z[1])",
    )
    GMT.grdimage(grd; xaxis = "af+l@[y@[", yaxis = "af+l@[z@[")
    GMT.plot!([sp.y sp.z sp.r * 14 * 2 / (last(nf.y) - first(nf.y))]; symbol = "c")
    GMT.colorbar!(;
        xaxis = "af",
        yaxis = "+l@[\\mathbf{Re}E_x@[",
        nolines = true,
        show = true,
        savefig = "Ex_re.png",
    )
end

# ╔═╡ 048f7096-6da0-445d-b690-69224bc67276
let
    spheres = param.spheres
    x = map(s -> s.x, spheres)
    y = map(s -> s.y, spheres)
    z = map(s -> s.z, spheres)
    r = map(s -> s.r, spheres)
    GMT.plot3d(
        [x y z fill(-20, length(x)) r 2r];
        symbol = "e",
        show = true,
        savefig = "setting.png",
    )
end

# ╔═╡ Cell order:
# ╠═93737264-8d97-11ec-2d60-336b353c2e54
# ╠═9e0a299f-9c69-4cb1-8de7-62116eb4e404
# ╠═4f69842a-3658-4345-a6ce-8dcd0c20c564
# ╠═a826a137-7c2c-4c70-9217-4797ca95b593
# ╠═048f7096-6da0-445d-b690-69224bc67276
