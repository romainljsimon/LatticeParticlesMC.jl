struct EXYZ <: Arianna.Format
    extension::String
    function EXYZ()
        return new(".exyz")
    end
end

function compute_box_str(box, ::EXYZ)
    if length(box) == 2
        return "$(box[1]) 0.0 0.0 0.0 $(box[2]) 0.0 0.0 0.0 0.0"
    elseif length(box) == 3
        return "$(box[1]) 0.0 0.0 0.0 $(box[2]) 0.0 0.0 0.0 $(box[3])"
    else
        throw(ArgumentError("Box vector must have 2 or 3 elements."))
    end
end
function write_header(io, system::Particles, t, format::EXYZ, digits::Integer)
    println(io, length(system))
    box_str = compute_box_str(system.box, format)
    println(io, "Lattice=\"$box_str\" Properties=$(get_system_column(system, format)):species:S:1:pos:R:$(system.d) Time=$t")
    return nothing
end
function get_system_column(::Molecules, ::EXYZ)
    return "molecule:I:1"
end
function formatted_string(num::Real, digits::Integer)
    fmtstr = "%." * string(digits) * "f"
    fmt = Printf.Format(fmtstr)
    return Printf.format(fmt, num)
end

function write_position(io, position, digits::Int)
    for position_i in position 
        formatted_position_i = formatted_string(position_i, digits)
        print(io, " ")
        print(io, formatted_position_i)
    end
    println(io)
    return nothing
end

function Arianna.store_trajectory(io, system::Molecules, t, format::Arianna.Format; digits::Integer=6)
    write_header(io, system, t, format, digits)
    for (molecule, species, position) in zip(system.molecule, system.species, system.position)
        print(io, "$molecule $species")
        write_position(io, position, digits)
    end
    return nothing
end
