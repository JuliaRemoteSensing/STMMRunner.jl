@inline function rand_point()
    return 2.0 .* rand(3) .- 1.0
end

function sample_in_sphere(R::Number, r::Number, n::Integer, max_retry_times::Integer = 100)
    r > R && throw(ArgumentError("particle radius must be smaller than region radius!"))

    points = Vector{Float64}[]
    r′ = r / R
    for _ = 1:n
        retry_times = 0
        new_point = rand_point()
        while norm(new_point) + r′ > 1.0
            new_point = rand_point()
        end
        min_dist = minimum(norm(point - new_point) for point in points; init = Inf)
        while min_dist < 2r′ && retry_times < max_retry_times
            retry_times += 1
            new_point = rand_point()
            while norm(new_point) + r′ > 1.0
                new_point = rand_point()
            end
            min_dist = minimum(norm(point - new_point) for point in points; init = Inf)
        end
        if min_dist >= 2r′
            push!(points, new_point)
        end
    end

    length(points) < n &&
        @warn "Failed to generate enough points. $(length(points))/$n generated."

    return [point ./ r′ for point in points]
end

export sample_in_sphere
