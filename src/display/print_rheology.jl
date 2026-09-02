# Pretty printing for CompositeRheologies
import Base.show

# returns a matrix with strings in the right order
function print_rheology_matrix(v::Tuple, el_num = nothing, digits = 1)

    n = 40
    A = Matrix{String}(undef, n, n)

    i, j, i_max = 1, 0, 1
    for entry in eachindex(v)
        out = print_rheology_matrix(v[entry], el_num[entry], digits)
        si = size(out)
        if prod(si) == 1
            j = j + 1
            A[i, j] = out[1]
        elseif (length(si) == 1) && prod(si) != 1
            j = j + 1
            A[i:(i + si[1] - 1), j] = out
            i_max = max(i_max, si[1])
        else
            A[i:(i + si[1] - 1), j:(j + si[2] - 1)] = out
            j = j + si[2] - 1
        end
    end
    for i in eachindex(A)
        if !isassigned(A, i)
            A[i] = ""
        end
    end
    for i in eachindex(A)
        if A[i] == ""
            A[i] = print_rheology_matrix("", digits)[1]
        end
    end

    A = A[1:i_max, 1:j]

    return A
end


function print_rheology_matrix(v::ParallelModel, el_num0 = nothing, digits = 1)
    n = 40
    A = Matrix{String}(undef, n, n)
    elements = superflatten((v.leafs, v.branches))
    if isnothing(el_num0)
        el_num0 = global_eltype_numbering(v)
        if maximum(superflatten(el_num0)) > 9
            digits = 2
        end
    end
    el_num = (el_num0[1]..., el_num0[2]...)
    i, j = 1, 1
    i_vec = Int[]
    for entry in eachindex(elements)
        out = print_rheology_matrix(elements[entry], el_num[entry], digits)
        si = size(out)
        if prod(si) == 1
            push!(i_vec, i)
            A[i, j] = out[1]
            i = i + 1
        elseif (length(si) == 1) && prod(si) != 1
            push!(i_vec, i)
            A[i:(i + si[1] - 1), j] = out
            i = i + si[1]
        else
            push!(i_vec, i)
            A[i:(i + si[1] - 1), j:(j + si[2] - 1)] = out
            i, j = i + si[1], j + si[2] - 1
        end
    end

    # Fill every index (with empty strings)
    for i in eachindex(A)
        if !isassigned(A, i)
            A[i] = ""
        end
    end

    # Extract the relevant part of A
    A = A[1:i, 1:j]

    # Center the strings & put brackets around it
    B = create_string_vec(A)

    #nel = maximum(textwidth.(remove_colors_string.(B)))
    #nel = maximum(textwidth.(B))
    nel = maximum(length_str_no_colors.(B))
    for i in 1:length(B)
        if any(in.(i_vec, i))
            #str_local = remove_colors_string(B[i])
            str_local = B[i]

            str_local = cpad(str_local, nel, "-")
            #str_local = B[i][1:5]*str_local*B[i][end-4:end]
            B[i] = "|" * str_local * "|"

        else
            B[i] = cpad(B[i], nel, " ")
            B[i] = "|" * B[i] * "|"
        end
    end

    return B
end

length_str_no_colors(str) = textwidth(remove_colors_string.(str)) + 1 * 6 + 2

function remove_colors_string(str::String)
    str = replace(str, r"\e\[[0-9;]*m" => "")
    return str
end

function create_string_vec(A)
    B = String[]
    for i in 1:size(A, 1)
        str1 = join(A[i, :])
        if length(str1) > 0
            push!(B, str1)
        end
    end

    return B
end

function print_rheology_matrix(v::SeriesModel, el_num0 = nothing, digits = 1)
    n = 40
    A = Matrix{String}(undef, n, n)
    elements = superflatten((v.leafs, v.branches))
    if isnothing(el_num0)
        el_num0 = global_eltype_numbering(v)
        if maximum(superflatten(el_num0)) > 9
            digits = 2
        end
    end
    el_num = (el_num0[1]..., el_num0[2]...)
    i, j, i_max = 1, 0, 1
    for entry in eachindex(elements)
        out = print_rheology_matrix(elements[entry], el_num[entry], digits)
        si = size(out)
        if prod(si) == 1
            j = j + 1
            A[i, j] = out[1]
        elseif (length(si) == 1) && prod(si) != 1
            j = j + 1
            A[i:(i + si[1] - 1), j] = out
            i_max = max(i_max, si[1])
        else
            j = j + 1
            A[i:(i + si[1] - 1), j:j] = out
            i = i + si[1]

            i_max = max(i_max, si[1])
        end
    end
    # Fill every index (with empty strings)
    for i in eachindex(A)
        if !isassigned(A, i)
            A[i] = ""
        end
    end

    for i in eachindex(A)
        if A[i] == ""
            A[i] = print_rheology_matrix("", digits)[1]
        end
    end

    A = A[1:i_max, 1:j]             # Extract the relevant part of A
    B = create_string_vec(A)        # Create strings

    return B
end

# Print the individual rheological elements in the REPL
# Note: would probably be good to define AbstractViscosity, AbstractElasticity in addition to AbstractPlasticity
# colors:
#           \e[34m - blue (for compressible elements)
#           \e[39m - default

print_rheology_matrix(v::String, digits = 1) = ["\e[39m  $(emptysuperscript(digits))       \e[39m"]
print_rheology_matrix(v::AbstractViscosity, n = 1, digits = 1) = ["\e[39m--⟦▪̲̅▫̲̅▫̲̅▫̲̅$(superscript(n, digits))--\e[39m"]
print_rheology_matrix(v::AbstractRheology, n = nothing, digits = 1) = ["\e[39m--?????$(superscript(n, digits))--\e[39m"]
function print_rheology_matrix(v::AbstractElasticity, n = nothing, digits = 1)
    if _isvolumetric(v)
        return ["\e[34m--/\\/\\/$(superscript(n, digits))--\e[39m"]
    else
        return ["\e[39m--/\\/\\/$(superscript(n, digits))--\e[39m"]
    end
end
print_rheology_matrix(v::AbstractPlasticity, n = nothing, digits = 1) = ["\e[39m--▬▬▬__$(superscript(n, digits))--\e[39m"]

# Center strings
cpad(s, n::Integer, p = " ") = rpad(lpad(s, div(n + textwidth(s), 2), p), n, p)

"""
Creates a superscript string for the given integer.
"""
function superscript(n, digits = 1)
    str = "$(n[1])"
    if n[1] < 10 && digits == 2
        str = "0" * str
    end
    str = replace(str, "0" => "⁰", "1" => "¹", "2" => "²", "3" => "³", "4" => "⁴", "5" => "⁵", "6" => "⁶", "7" => "⁷", "8" => "⁸", "9" => "⁹")

    return str
end

superscript(n::Nothing, digits = 1) = ""
function emptysuperscript(digits = 1)
    if digits == 1
        str = " "
    elseif digits == 2
        str = "  "
    end
    return str
end


# Print the individual rheological elements in the REPL
function Base.show(io::IO, c::SeriesModel)
    str = print_rheology_matrix(c)
    println.(io, str)
    return nothing
end

function Base.show(io::IO, c::ParallelModel)
    str = print_rheology_matrix(c)
    println.(io, str)
    return nothing
end
