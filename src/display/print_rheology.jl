# Pretty printing for CompositeRheologies

# Upper bound on the rows and columns a rendered block can occupy: each rheology
# element contributes at most one row and one column, and nesting only
# redistributes them between the two. Sizing the buffer from this rather than
# from a fixed constant means an oversized composite raises a BoundsError
# instead of being silently trimmed.
_n_rendered(::AbstractRheology) = 1
_n_rendered(c::AbstractCompositeModel) = _n_rendered(c.leafs) + _n_rendered(c.branches)
_n_rendered(t::Tuple) = sum(_n_rendered, t; init = 0)

# Elements of a composite in the order they are drawn, with their per-type
# numbering. Two digits are needed as soon as any counter passes 9.
function _rendered_elements(v::AbstractCompositeModel, el_num0, digits)
    elements = superflatten((v.leafs, v.branches))
    if isnothing(el_num0)
        el_num0 = global_eltype_numbering(v)
        maximum(superflatten(el_num0)) > 9 && (digits = 2)
    end
    return elements, (el_num0[1]..., el_num0[2]...), digits
end

"""
    _place_blocks(v, el_num0, digits, grow_rows)

Render each element of `v` into its own block of text and lay the blocks out in
a matrix, one column per element for a series composite (`grow_rows = false`,
elements drawn left to right) and one row-span per element for a parallel one
(`grow_rows = true`, elements stacked). Unassigned cells become empty strings.

Returns the trimmed matrix, the row each block starts on (parallel only, for
bracketing), and the digit width the numbering settled on.
"""
function _place_blocks(v::AbstractCompositeModel, el_num0, digits, grow_rows::Bool)
    elements, el_num, digits = _rendered_elements(v, el_num0, digits)
    n = _n_rendered(v) + 1
    A = Matrix{String}(undef, n, n)

    i, j, height = 1, (grow_rows ? 1 : 0), 1
    starts = Int[]
    for entry in eachindex(elements)
        block = print_rheology_matrix(elements[entry], el_num[entry], digits)
        rows = length(block)
        if grow_rows
            push!(starts, i)
            A[i:(i + rows - 1), j] = block
            i += rows
        else
            j += 1
            A[i:(i + rows - 1), j] = block
            height = max(height, rows)
        end
    end

    for k in eachindex(A)
        isassigned(A, k) || (A[k] = "")
    end

    return A[1:(grow_rows ? i : height), 1:j], starts, digits
end

function print_rheology_matrix(v::ParallelModel, el_num0 = nothing, digits = 1)
    A, starts, _ = _place_blocks(v, el_num0, digits, true)
    B = create_string_vec(A)

    # Bracket each row. Rows a block starts on are filled with dashes, so the
    # branch reads as connected; the rest are filled with spaces.
    width = maximum(length_str_no_colors.(B))
    for i in eachindex(B)
        fill = i in starts ? "-" : " "
        B[i] = "|" * cpad(B[i], width, fill) * "|"
    end

    return B
end

# Columns of dash- or space-fill that `print_rheology_matrix(::ParallelModel)`
# lays around each row of a bracketed block. It sets the width of that fill and
# nothing else.
const ROW_FILL_WIDTH = 8

# Visible width of `str`: `textwidth` counts ANSI colour escapes, which occupy no
# columns, so they are stripped first.
length_str_no_colors(str) = textwidth(remove_colors_string(str)) + ROW_FILL_WIDTH

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
    A, _, digits = _place_blocks(v, el_num0, digits, false)

    # Blank cells become an element-sized run of spaces, so columns stay aligned
    # where one branch is taller than its neighbours.
    blank = print_rheology_matrix("", digits)[1]
    for k in eachindex(A)
        A[k] == "" && (A[k] = blank)
    end

    return create_string_vec(A)
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
function Base.show(io::IO, c::Union{SeriesModel, ParallelModel})
    println.(io, print_rheology_matrix(c))
    return nothing
end
