module LatticeParticlesMC

using Arianna
using StaticArrays
using Statistics
using Printf

export Particles
abstract type Particles <: AriannaSystem end
# Write your package code here.
include("molecules.jl")
include("utils.jl")
include("moves.jl")
include("io.jl")
get_position(system::Particles, i::Int) = @inbounds system.position[i]
get_species(system::Particles, i::Int) = @inbounds system.species[i]
get_model(system::Particles, i::Int, j::Int) = @inbounds system.model_matrix[get_species(system, i), get_species(system, j)]
get_box(system::Particles) = system.box
get_neighbour_list(system::Particles) = system.neighbour_list
Base.length(system::Particles) = system.N
Base.eachindex(system::Particles) = Base.OneTo(length(system))

export System, Molecules
export SimpleUniform
export Displacement, Corner
export callback_energy
export EXYZ
end
