function Arianna.unnormalised_log_target_density(e, system::Particles)
    return  - e / system.temperature
end

function Arianna.delta_log_target_density(e1, e2, system::Particles)
    return -(e2 - e1) ./ system.temperature
end

fold_back(x, box) = x .- fld.(x, box) .* box


function vector_1D(c1, c2, side_length)
    dx = c1 - c2
    return dx - round(dx / side_length) * side_length
end

@inline function vector(c1::SVector{N,Int}, c2::SVector{N,Int}, box::SVector{N,Int}) where {N}
    @inbounds return SVector{N,Int}(ntuple(i -> vector_1D(c1[i], c2[i], box[i]), Val(N)))
end

@inline function nearest_image_distance_squared(xi::SVector{N,Int}, xj::SVector{N,Int},
                                                box::SVector{N, Int}) where {N}
    @inbounds dx = vector(xi, xj, box)
    return sum(abs2, dx)
end
