module Sequences
using Combinatorics

"""
    find_index(seq1::Vector{<:Any}, seq2::Vector{<:Any})

Find the index of each of the members of seq1 in seq2
"""
function find_index(seq1::Vector{<:Any}, seq2::Vector{<:Any})
    idx = fill(0, length(seq1))
    for (ii,s1) in enumerate(seq1)
        jj = findfirst(seq2.==s1)
        if jj !== nothing
            idx[ii] = jj
        end
    end
    idx
end

function count_occurrences(sub::Vector{T}, sequences::Vector{Vector{T}}) where T
    count = 0
    for seq in sequences
        if is_subsequence(sub, seq)
            count += 1
        end
    end
    return count
end

function is_subsequence(sub::Vector{T}, seq::Vector{T}) where T
    sub_idx = 1
    for val in seq
        if sub_idx <= length(sub) && val == sub[sub_idx]
            sub_idx += 1
        end
    end
    return sub_idx > length(sub)
end

function longest_frequent_subsequence(sequences::Vector{Vector{T}}, min_sequences=1) where T
    # Find the shortest sequence to minimize the search space
    shortest_seq_idx = argmin(length.(sequences))
    shortest_seq = sequences[shortest_seq_idx]
    n = length(shortest_seq)
    best_subseq = T[]
    best_count = 0

    # Generate all possible subsequences of the shortest sequence, from longest to shortest
    for len in n:-1:1
        for indices in combinations(1:n, len)
            candidate = shortest_seq[indices]
            count = count_occurrences(candidate, sequences)
            if count >= min_sequences && (length(candidate) > length(best_subseq) || (length(candidate) == length(best_subseq) && count > best_count))
                best_subseq = candidate
                best_count = count
            end
        end
    end
    return best_subseq, best_count
end

function lcs_multiple(sequences::Vector{Vector{T}}) where T
    # Find the shortest sequence to minimize the search space
    shortest_seq_idx = argmin(length.(sequences))
    shortest_seq = sequences[shortest_seq_idx]
    n = length(shortest_seq)
    best_lcs = T[]

    # Generate all possible subsequences of the shortest sequence, from longest to shortest
    for len in n:-1:1
        for indices in combinations(1:n, len)
            candidate = shortest_seq[indices]
            # Check if candidate is a subsequence of all sequences
            all_good = true
            for seq in sequences
                if !is_subsequence(candidate, seq)
                    all_good = false
                    break
                end
            end
            if all_good
                return candidate
            end
        end
    end
    return best_lcs
end

mutable struct State
    len::Int          # Length of the longest substring in this state
    link::Int         # Suffix link
    next::Dict{Int,Int}  # Transitions (character -> state)
end

function build_suffix_automaton(s::Vector{Int})
    sa = [State(0, -1, Dict{Int,Int}())]
    last = 1
    for c in s
        p = last
        curr = length(sa) + 1
        push!(sa, State(sa[p].len + 1, -1, Dict{Int,Int}()))
        while p != -1 && !haskey(sa[p].next, c)
            sa[p].next[c] = curr
            p = sa[p].link
        end
        if p == -1
            sa[curr].link = 1
        else
            q = sa[p].next[c]
            if sa[p].len + 1 == sa[q].len
                sa[curr].link = q
            else
                clone = length(sa) + 1
                push!(sa, State(sa[p].len + 1, sa[q].link, copy(sa[q].next)))
                while p != -1 && sa[p].next[c] == q
                    sa[p].next[c] = clone
                    p = sa[p].link
                end
                sa[q].link = clone
                sa[curr].link = clone
            end
        end
        last = curr
    end
    return sa
end

function longest_common_contiguous_subsequence(s1::Vector{Int}, s2::Vector{Int})
    sa = build_suffix_automaton(s1)
    max_len = 0
    max_pos = 0
    v = 1
    l = 0
    for (i, c) in enumerate(s2)
        while v != -1 && !haskey(sa[v].next, c)
            v = sa[v].link
            l = (v == -1) ? 0 : sa[v].len
        end
        if v == -1
            v = 1
            l = 0
            continue
        end
        v = sa[v].next[c]
        l += 1
        if l > max_len
            max_len = l
            max_pos = i - l + 1
        end
    end
    return s2[max_pos:max_pos+max_len-1]
end

function lccs_multiple(sequences::Vector{Vector{Int}})
    if isempty(sequences)
        return Int[]
    end
    current_lccs = sequences[1]
    for seq in sequences[2:end]
        current_lccs = longest_common_contiguous_subsequence(current_lccs, seq)
        isempty(current_lccs) && break
    end
    return current_lccs
end

end