using Distributions

@inline function rand_point()
    return 2.0 .* rand(3) .- 1.0
end

function sample_in_sphere(rfunc, R::Number, n::Integer, max_retry_times::Integer=100)
    points = Vector{Float64}[]
    radius = Float64[]
    for _ = 1:n
        retry_times = 0
        r = rfunc()
        while r > R
            r = rfunc()
        end
        new_point = rand_point() .* (R - r)

        min_dist = minimum(norm(point - new_point) - r - r₀ for (r₀, point) in zip(radius, points); init=Inf)
        while min_dist < 0 && retry_times < max_retry_times
            retry_times += 1
            new_point = rand_point() .* (R - r)
            while norm(new_point) + r > R
                r = rfunc()
                while r > R
                    r = rfunc()
                end
                new_point = rand_point() .* (R - r)
            end
            min_dist = minimum(norm(point - new_point) - r - r₀ for (r₀, point) in zip(radius, points); init=Inf)
        end
        if min_dist >= 0
            push!(points, new_point)
            push!(radius, r)
        end
    end

    length(points) < n &&
        @warn "Failed to generate enough points. $(length(points))/$n generated."

    return [vcat(point, r) for (r, point) in zip(radius, points)]
end

sample_in_sphere(R::Number, r::Number, n::Integer, max_retry_times::Integer=100) = sample_in_sphere(() -> r, R, n, max_retry_times)

sample_in_sphere_lognormal(R::Number, r::Number, σ::Number, n::Integer, max_retry_times::Integer=100) = sample_in_sphere(() -> rand(LogNormal(0, σ)) * r, R, n, max_retry_times)

export sample_in_sphere
