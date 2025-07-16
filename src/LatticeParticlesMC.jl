module LatticeParticlesMC

using Arianna, StaticArrays, Statistics, Printf

export Particles
abstract type Particles <: AriannaSystem end


include("utils.jl")

include("molecules.jl")
include("moves.jl")
include("io.jl")

get_position(system::Particles, i::Int) = @inbounds system.position[i]
get_species(system::Particles, i::Int) = @inbounds system.species[i]
get_model(system::Particles, i::Int, j::Int) = @inbounds system.model_matrix[get_species(system, i), get_species(system, j)]
get_box(system::Particles) = system.box
Base.length(system::Particles) = system.N
Base.eachindex(system::Particles) = Base.OneTo(length(system))

function compute_energy_particle(system::Particles, ids::AbstractVector)
    return map(i -> compute_energy_particle(system, i), ids)
end

export callback_energy
export Molecules
export Displacement, Corner
export System
export SimpleUniform
export EXYZ
end
