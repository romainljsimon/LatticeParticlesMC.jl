struct Molecules{D,  VS<:AbstractVector,  T<:AbstractFloat, SM<:AbstractArray} <: Particles
    position::Vector{SVector{D,Int}}
    species::VS
    molecule::VS
    molecule_species::VS
    start_mol::Vector{Int}
    length_mol::Vector{Int}
    density::T
    temperature::T
    energy::Vector{T}
    model_matrix::SM
    d::Int
    N::Int
    Nmol::Int
    box::SVector{D,Int}
    bonds::Vector{Vector{Int}}
end

function System(position, species, molecule, box, temperature::T, model_matrix, bonds; molecule_species=nothing) where {T<:AbstractFloat}
    @assert length(position) == length(species)
    N = length(position)
    Nmol = length(unique(molecule))
    start_mol, length_mol = get_first_and_counts(molecule)
    molecule_species = something(molecule_species, ones(Int, N))
    d = length(Array(position)[1])
    #box = @SVector fill(T((N / density)^(1 / d)), d)
    density = N / prod(box)
    energy = zeros(T, 1)
    system = Molecules(position, species, molecule, molecule_species,  start_mol, length_mol, density, temperature, energy, model_matrix, d, N, Nmol, box, bonds)
    local_energy = [compute_energy_particle(system, i) for i in eachindex(position)]
    energy = sum(local_energy) / 2
    if isinf(energy) || isnan(energy)
        error("Initial configuration has infinite or NaN energy.")
    end
    system.energy[1] = energy
    return system
end

get_start_end_mol(system::Molecules, i::Int) = system.start_mol[i], system.start_mol[i] + system.length_mol[i] - 1

function get_first_and_counts(vec::Vector{Int})
    firsts = Int[]
    counts = Int[]
    
    # Handle empty vector case
    isempty(vec) && return firsts, counts
    
    # Initialize with first element
    current = vec[1]
    push!(firsts, 1)
    count = 1
    
    # Scan through vector
    @inbounds for i in 2:length(vec)
        if vec[i] != current
            push!(counts, count)
            push!(firsts, i)
            current = vec[i]
            count = 1
        else
            count += 1
        end
    end
    # Add last count
    push!(counts, count)
    
    return firsts, counts
end

function compute_energy(system::Particles)
    energy = zero(typeof(system.temperature))
    for i in 1:system.N
        energy += compute_energy_particle(system, i)
    end
    return energy / 2
end

function compute_energy_outer_molecule(system, molecule)
    energy = zero(typeof(system.temperature))
    for molecule_j in 1:system.Nmol
        molecule_j == molecule && continue
        energy += compute_energy_molecule_molecule(system, molecule, molecule_j)
    end
    return energy
end


function compute_energy_molecule_molecule(system::Molecules, molecule_i, molecule_j)
    start_mol_i, end_mol_i = get_start_end_mol(system, molecule_i)
    start_mol_j, end_mol_j = get_start_end_mol(system, molecule_j)
    energy = zero(typeof(system.temperature))
    for i in start_mol_i:end_mol_i
        for j in start_mol_j:end_mol_j
            position_i, position_j = system.position[i], system.position[j]
            interaction_ij = system.model_matrix[i, j]
            energy += compute_energy_ij(system, position_i, position_j, interaction_ij)
        end
    end
    return energy
end

function check_compute_energy_ij(system::Molecules, i, j, position_i, bonds_i)
    # Early return using && for short-circuit evaluation
    i == j && return zero(typeof(system.temperature))
    j ∈ bonds_i && return zero(typeof(system.temperature))
    # If i != j, compute energy directly
    position_j = get_position(system, j)
    interaction_ij = get_model(system, i, j)
    return compute_energy_ij(system, position_i, position_j, interaction_ij)
end

function compute_energy_ij(system::Molecules, position_i, position_j, interaction_ij)
    r2 = nearest_image_distance_squared(position_i, position_j, get_box(system))
    r2 == 0 && return Inf
    cutoff2_val = 1
    return r2 > cutoff2_val ? zero(typeof(system.temperature)) : interaction_ij
end

function compute_energy_particle(system::Molecules, i)
    energy = zero(typeof(system.temperature))
    position_i = system.position[i]
    bonds_i = system.bonds[i]
    for j in eachindex(system)
        energy += check_compute_energy_ij(system, i, j, position_i, bonds_i)
    end
    return energy
end

function callback_energy(simulation)
    return mean(system.energy[1] / length(system) for system in simulation.chains)
end
