using LatticeParticlesMC
using Arianna
using StaticArrays

using Distributions
using Random

seed = 42
rng = Xoshiro(seed)
N = 25
M = 1
d = 3
temperature = 1.0
density = 1.2

function create_bond_matrix(N::Int)
    # Create a vector to store the SVectors, each containing a pair of integers
    matrix = Vector{Vector{Int}}()
    # Populate the matrix with pairs according to the specified pattern
    push!(matrix, [2])
    for i in 2:N-1
        push!(matrix, [i-1, i + 1])
    end
    push!(matrix, [N-1])
    return matrix
end


bonds = create_bond_matrix(N)
position = [SVector{3, Int}(i, N÷2, N÷2) for i in 1:N]
species = collect(1:N)
molecule = [1 for i in 1:N]
interaction_matrix = randn(N, N)
interaction_matrix = 0.5 * (interaction_matrix + interaction_matrix')
box = SVector{3, Int}(N+2, N+2, N+2)
system = System(position, species, molecule, box, temperature, interaction_matrix, bonds)
chains=[system]
displacement_policy = SimpleUniform()
pool = (
    Move(Displacement(0, zero(chains[1].box), 0.0), displacement_policy, Vector{Float64}(), 0.1),
    Move(Corner(0, zero(chains[1].box), 0.0), displacement_policy, Vector{Float64}(), 0.9),
)
## Define the simulation struct
steps = 1000
burn = 0
block = [0, 10]
sampletimes = build_schedule(steps, burn, block)
callbacks = (callback_energy, callback_acceptance)

path = "data/test/particles/Molecules/T$temperature/N$N/M$M/seed$seed"
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=N),
    (algorithm=StoreCallbacks, callbacks=(callback_energy, callback_acceptance), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=EXYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10))
    )
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)