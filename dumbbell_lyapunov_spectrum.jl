#;cd ~

#;cd /home/rafael/UNAM/Sanders/Mancuerna/8-9-Agosto_y_Septiembre/Dumbbell-Julia

include("dumbbell_type.jl")
include("Gram_Schmidt.jl")
include("dumbbell_lyapunov_functions.jl")


using PyCall
pygui(:gtk)

using PyPlot

outfile = "exponents.dat"
f1 = open(outfile,"w")

rcm = [1., 1.]; l0=0.2; m0=1; vcm = [2.1, 0.6]; θ=π/6; ω=2

vert_walls = [-3., 3.]; horiz_walls =[-2., 2.]

db = dumbbell(rcm, θ, vcm, ω, m0, l0, 0)

#reset_dummbell(db)

randvecs = ([rand(6) for i in 1:6])
dispvecs, normdispvecs = Gram_Schmidt(randvecs)

N=1000

col_count = 0

λ_col = zeros(Float64,N,6)
λ_time = zeros(Float64,N,6)
χs = zeros(Float64,N,6)

col_times = Float64[]; times = Float64[]; t_T = 0

check=true

#tic()
for k in 1:N
    col_count += 1
    
    Γ0,  Γf, wall, part, tc = collision(db, vert_walls, horiz_walls)
    t_T += tc
    push!(times, t_T)
    push!(col_times, tc)
    
    dispvecs0 = Array{Float64}[]
    for v in normdispvecs
        vf = evol_Γ(v, tc, l0, m0)
        push!(dispvecs0, [vf])
    end
    
    dispvecs1 = Array{Float64}[]
    for dΓ in dispvecs0
        dΓp = displacement_vector_collision_map(Γ0, dΓ, part, wall, l0, m0)
        push!(dispvecs1, [dΓp])
    end
    
    orthovecs, normdispvecs = Gram_Schmidt(dispvecs1)
    
    for (j,vv) in enumerate(orthovecs)
        if k>1
            χs[k, j] = χs[k-1, j] + log(norm(vv))
        else
            χs[1,j] = log(norm(vv))
        end
        
        λ_col[k,j] = χs[k,j]/(k+1)
        λ_time[k,j] = χs[k,j]/t_T
    end
    
    println(f1, string(t_T)*"\t"*string(tc))
    for j in 1:6
        print(f1,"\t"*string(λ_col[k,j])*"\t"*string(λ_time[k,j]))
    end
    print(f1,"\t"*string(db.collision_counter))
    
    x1, y1 = particle_position(db, 1)[1:2]
    x2, y2 = particle_position(db, 2)[1:2]
    
    if k%1000==0 && check
        suma1 = 0.; suma2 = 0.
        for j in 1:6
            suma1 += λ_time[k,j]
            suma2 += λ_col[k,j]
        end
        
        println("Suma de exponentes; \t"*string(db.collision_counter)*" colisiones")
        println(suma1, "\t", suma2)
        println("--------\nDumbbell position")
        println(x1,"\t", y1,"\t", x2,"\t", y2)
        #toc()
        println("*********\n")
        #tic()
    end
    
end
    
close(f1)

println("Total time:\t", t_T, "\n********")

println("Lyapunov Spectrum")
for j in 1:6
    println(λ_col[end,j], "\t", λ_time[end,j], "\t", λ_col[end,j]/λ_time[end,j])
end
println("Tiempo promedio \t",t_T/N)


colors = ["cyan", "blue", "red", "green", "purple", "yellow"]
figure(figsize=(10,10))

for k in 1:6
    plot(λ_col[:,k], color=colors[k])
    plot(λ_time[:,k], color=colors[k], ls="--")
end

show()