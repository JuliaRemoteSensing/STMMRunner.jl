### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ a343426c-952d-11ec-0693-3dbb707e1e76
begin
    import Pkg
    Pkg.activate("..")

    import GMT

    using DelimitedFiles
    using STMMRunner
    using STMMRunner.MSTM.V4
end

# ╔═╡ 1b80ded1-6068-4370-9aa7-3f04ab013b11
begin
    spheres = [
		Sphere(
			r = 0.7 + 0.5 * rand(),
			x = i,
			y = j,
			z = k,
			m = 1.5 + 1.0im * rand()
		) for i = -6.0:3.0:6.0 for j = -6.0:3.0:6.0 for k = -6.0:3.0:6.0
	]

   param = STMMConfig(
		number_processors = Sys.CPU_THREADS ÷ 2,
		sm_number_processors = Sys.CPU_THREADS ÷ 2,
		run_print_file = "",
		working_directory = "",
		orientation = FixedOrientation,
		α = 70.0,
		calculate_near_field = true,
		near_field_x_min = 0.0,
		near_field_y_min = -10.0,
	    near_field_z_min = -10.0,
	    near_field_x_max = 0.0,
	    near_field_y_max = 10.0,
	    near_field_z_max = 10.0,
	    near_field_step_size = 0.1,
		spheres = spheres,
	)
end

# ╔═╡ ffcb61ce-843f-4108-9942-bc00ae3ca36a
results = run_mstm(param; keep = false)

# ╔═╡ 7ef49694-b211-4f66-9c29-6b9704c775fc
let
    nf = results.near_field.field
    y = nf.y |> Set |> collect |> sort
    z = nf.z |> Set |> collect |> sort
    sp = results.near_field.spheres

    grd = GMT.xyz2grd(
        [nf.y nf.z real.(nf.Ex₌)],
        region = (first(y), last(y), first(z), last(z)),
        inc = "$(y[2]-y[1])/$(z[2]-z[1])",
    )
    GMT.grdimage(grd; xaxis = "af+l@[y@[", yaxis = "af+l@[z@[")
    GMT.plot!([sp.y sp.z sp.r * 14 * 2 / (last(nf.y) - first(nf.y))]; symbol = "c")
    GMT.colorbar!(;
        xaxis = "af",
        yaxis = "+l@[\\mathbf{Re}E_{{\\parallel}x}@[",
        nolines = true,
        show = true,
        savefig = "E∥x_re_v4.png",
    )
end

# ╔═╡ 224530e6-1fee-44b0-a7d8-be879224b361
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
        savefig = "setting_v4.png",
    )
end

# ╔═╡ Cell order:
# ╠═a343426c-952d-11ec-0693-3dbb707e1e76
# ╠═1b80ded1-6068-4370-9aa7-3f04ab013b11
# ╠═ffcb61ce-843f-4108-9942-bc00ae3ca36a
# ╠═7ef49694-b211-4f66-9c29-6b9704c775fc
# ╠═224530e6-1fee-44b0-a7d8-be879224b361
