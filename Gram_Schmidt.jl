
#Function for finding the projection of the 'u' vector along the 'v' one.

function proj(u::Vector,v::Vector)
    l = length(u)
    if length(u[u.!=0])==0
        return zeros(Float64,l)
    else
        return (u⋅v/(norm(u)^2))*u
    end
end


# Function for orthogonalize and orthonormalize the set of vectors 'vecs' using the Grand-Schmidt procedure
# The functions returns two sets of vectors, the orthogonal ones and the orthoNORMAL ones.
function Gram_Schmidt(vecs)
    N = length(vecs) #number of vectors
    dim = length(vecs[1]) #dimension of vectors
    v0 = vecs[1]
    ovecs = Array{Float64}[]
    nvecs = Array{Float64}[]

    for k in 1:N
        projections = zeros(Float64,dim)
        for u in ovecs
            projections += proj(u, vecs[k])
        end
        push!(ovecs, [vecs[k]-projections])
        push!(nvecs, [ovecs[k]/norm(ovecs[k])])
    end
    ovecs, nvecs
end