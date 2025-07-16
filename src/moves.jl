
function Arianna.perform_action!(system::Particles, action::Action)
    e₁, e₂ = perform_action!(system, action)
    if isinf(e₁) || isinf(e₂)
        action.δe = zero(typeof(system.energy[1]))
    else
        action.δe = e₂ - e₁
        system.energy[1] += action.δe
    end
    return e₁, e₂
end

###############################################################################
# SIMPLE DISPLACEMENT

"""
    mutable struct Displacement{T<:AbstractArray} <: Action

A struct representing a displacement action, where particle i is moved by specified amounts `δ`.

# Fields
- `i::Int`: Indices of the particles or elements in `system` to be displaced.
- `δ::T`: Displacement values for each corresponding index in `is`.
"""
mutable struct Displacement{T<:AbstractArray, F<:AbstractFloat} <: Action
    i::Int
    δ::T
    δe::F
end

function update_position!(system::Particles, action::Displacement)
    @inbounds system.position[action.i] = system.position[action.i] + action.δ
end

function perform_action!(system::Particles, action::Displacement)
    e₁ = compute_energy_particle(system, action.i)
    update_position!(system, action)
    e₂ = compute_energy_particle(system, action.i)
    action.δe = e₂ - e₁
    return e₁, e₂
end

function Arianna.invert_action!(action::Displacement, ::Particles)
    action.δ = -action.δ
end


struct SimpleUniform <: Policy end

function Arianna.log_proposal_density(::Displacement, ::SimpleUniform, parameters, system::Particles)
    return -1
end

function Arianna.sample_action!(action::Displacement, ::SimpleUniform, parameters, system::Particles, rng)
    # Select on molecule, then select one end atom and move it
    molecule = rand(rng, 1:system.Nmol)
    start_mol, end_mol = get_start_end_mol(system, molecule)

    action.i = rand(rng, (start_mol, end_mol))
    if action.i == start_mol
        direction_index = action.i + 1
    else
        direction_index = action.i - 1
    end
        direction = system.position[action.i] - system.position[direction_index]
    i = rand(1:system.d)              # choose axis
    s = rand((1, -1))          # choose direction
    v = zeros(Int, system.d)          # base vector
    v[i] = s

    action.δ = -direction + v
end



###############################################################################
# SIMPLE DISPLACEMENT

"""
    mutable struct Displacement{T<:AbstractArray} <: Action

A struct representing a displacement action, where particle i is moved by specified amounts `δ`.

# Fields
- `i::Int`: Indices of the particles or elements in `system` to be displaced.
- `δ::T`: Displacement values for each corresponding index in `is`.
"""
mutable struct Corner{T<:AbstractArray, F<:AbstractFloat} <: Action
    i::Int
    δ::T
    δe::F
end

function update_position!(system::Particles, action::Corner)
    @inbounds system.position[action.i] = system.position[action.i] + action.δ
end

function perform_action!(system::Particles, action::Corner)
    e₁ = compute_energy_particle(system, action.i)
    update_position!(system, action)
    e₂ = compute_energy_particle(system, action.i)
    action.δe = e₂ - e₁
    return e₁, e₂
end

function Arianna.invert_action!(action::Corner, ::Particles)
    action.δ = -action.δ
end


function Arianna.log_proposal_density(::Corner, ::SimpleUniform, parameters, system::Particles)
    return -1
end

function Arianna.sample_action!(action::Corner, ::SimpleUniform, parameters, system::Particles, rng)
    # Select on molecule, then select one end atom and move it
    molecule = rand(rng, 1:system.Nmol)
    if system.length_mol[molecule] < 2
        action.δ = 0.0
        return
    end
    start_mol, end_mol = get_start_end_mol(system, molecule)
    action.i = rand(rng, start_mol+1:end_mol-1)
    direction1 = system.position[action.i] - system.position[action.i-1]
    direction2 = system.position[action.i] - system.position[action.i+1]
    action.δ = -direction1 - direction2
end

