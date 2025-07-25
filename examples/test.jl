using LatticeParticlesMC
using Arianna
using Random
using StaticArrays
using DelimitedFiles

seed = 42
rng = Xoshiro(seed)
N1 = 25
N2 = 5
N = N1 + N2
M = 1
d = 3
temperature = 1.0
box = SVector{d, Int}(N+1, N+1, N+1)

function create_bond_matrix(N::Int; step::Int=0)
    # Create a vector to store the SVectors, each containing a pair of integers
    matrix = Vector{Vector{Int}}()
    # Populate the matrix with pairs according to the specified pattern
    push!(matrix, [step+2])
    for i in 2:N-1
        push!(matrix, [step+i-1, step+i + 1])
    end
    push!(matrix, [step+N-1])
    return matrix
end

#position1 = [SVector{d, Int}(i, N÷2, N÷2) for i in 1:N1]
#position2 = [SVector{d, Int}(i, 1, 1) for i in 1:N2]
#position = [position1; position2]

bonds1 = create_bond_matrix(N1)
bonds2 = create_bond_matrix(N2; step=25)
bonds = [bonds1; bonds2]

#molecule1 = [1 for i in 1:N1]
#molecule2 = [2 for i in 1:N2]
#molecule = [molecule1; molecule2]

#species = collect(1:N)


lines = readlines("examples/init.xyz")

# Skip the first 2 lines (header and lattice line)
data_lines = lines[3:end]

# Prepare empty vectors
molecule = Int[]
species = Int[]
position = Vector{SVector{3, Int}}()  # or SVector if using StaticArrays

# Parse each data line
for line in data_lines
    fields = split(line)
    push!(molecule, parse(Int, fields[1]))
    push!(species, parse(Int, fields[2]))
    pos = (parse(Int, fields[3]), parse(Int, fields[4]), parse(Int, fields[5]))
    push!(position, pos)
end

interaction_matrix = readdlm("examples/interaction.dat")


molecules = System(position, species, molecule, box, temperature, interaction_matrix, bonds)
chains = [molecules]

policy = SimpleUniform()
parameters = Vector{Float64}()
pool = (
    Move(Tail(0, zero(chains[1].box), 0.0), policy, parameters, 0.1),
    Move(Displacement(0, zero(chains[1].box), 0.0), policy, parameters, 0.1),
    Move(Corner(0, zero(chains[1].box), 0.0), policy, parameters, 0.7 ),
    Move(Crankshaft(0, zero(chains[1].box), 0.0), policy, parameters, 0.1 ),
)
## Define the simulation struct
steps = 100000
burn = 0
block = [0, 10000]
sampletimes = build_schedule(steps, burn, block)
callbacks = (callback_acceptance)

path = "data/test/particles/Molecules/T$temperature/N$N/M$M/seed$seed"
algorithm_list = (
    (algorithm=Metropolis, pool=pool, seed=seed, parallel=false, sweepstep=N),
    (algorithm=StoreCallbacks, callbacks=(callback_acceptance, callback_energy), scheduler=sampletimes),
    (algorithm=StoreTrajectories, scheduler=sampletimes, fmt=EXYZ()),
    (algorithm=PrintTimeSteps, scheduler=build_schedule(steps, burn, steps ÷ 10))
    )
simulation = Simulation(chains, algorithm_list, steps; path=path, verbose=true)
run!(simulation)