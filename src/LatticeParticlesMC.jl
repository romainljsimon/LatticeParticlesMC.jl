module LatticeParticlesMC

using Arianna
using StaticArrays
using Statistics
using Printf
using LinearAlgebra

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
Base.length(system::Particles) = system.N
Base.eachindex(system::Particles) = Base.OneTo(length(system))

export System, Molecules
export SimpleUniform
export Tail, Corner, Displacement, Crankshaft
export callback_energy
export nearest_image_distance_squared
export EXYZ
end
