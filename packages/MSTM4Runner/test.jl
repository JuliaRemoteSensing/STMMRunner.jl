using MSTM4Runner

sphere1 = Sphere(r = 1.0, x = 4.2, y = 4.3, z = 4.0, m = 1.5)
sphere2 = Sphere(r = 1.0, x = 0.0, y = 0.3, z = 1.0, m = 1.5)

spheres = [sphere1, sphere2]
layers = [Layer(), Layer(d = 10.0, m = 2.0), Layer(m = 1.2)]

base = STMMConfig(number_processors = Sys.CPU_THREADS ÷ 2,
                  scattering_map_model = 0,
                  α = 5.0,
                  β = 12.0,
                  spheres = spheres,
                  frame = IncidentFrame)

fixed_base = STMMConfig(base;
                        orientation = FixedOrientation,
                        calculate_near_field = true,
                        near_field_x_min = -5.0,
                        near_field_x_max = 5.0,
                        near_field_y_min = -5.0,
                        near_field_y_max = 5.0,
                        near_field_z_min = -5.0,
                        near_field_z_max = 5.0,
                        near_field_step_size = 1.0)

param = STMMConfig(fixed_base)
run_mstm(param; keep = false)
