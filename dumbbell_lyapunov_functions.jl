l0=0.2; m0=1.;

function coefficients(θ, part, wall, l=l0)
    if wall=="horizontal"
        if part==1
            α = π/2 - θ
        elseif part==2
            α = -(π/2 - θ)
        end
    elseif wall=="vertical"
        if part ==1
            α = -θ
        elseif part==2
            α = θ
        end
    end
    a = 2/(1+csc(α)^2) - 1
    b = 2*l*csc(α)/(1+csc(α)^2)
    c = b/l^2
    a, b, c
end


function collision_map(Γ, part, wall, l=l0)
    qq = Γ[1:3]
    x, y, θ = qq
    pp = Γ[4:6]
    px, py, L = pp

    a, b, c, = coefficients(θ, part, wall, l)

    if wall=="horizontal"
        n1=0; n2=1;
    elseif wall=="vertical"
        n1=1; n2=0;
    end

    B_mat =
    [n2^2+a*n1^2 n1*n2*(a-1) -c*n1;
    n1*n2*(a-1) n1^2+a*n2^2 -c*n2;
    -b*n1 -b*n2 -a]

    pp = B_mat*pp

    return [qq, pp]
end


function collision_map_derivative(Γ, vec, part, wall, l=l0)
    qq = Γ[1:3];   x, y, θ = qq
    pp = Γ[4:6];   px, py, L = pp
    vecqq = vec[1:3]
    vecpp = vec[4:end]

    if part==1
        derα=-1
        if wall=="horizontal"
            α = π/2 - θ
        elseif wall=="vertical"
            α = -θ
        end

    elseif part==2
        derα = 1
        if wall=="horizontal"
            α = θ-π/2
        elseif wall=="vertical"
            α = θ
        end
    end

    a, b, c = coefficients(θ, part, wall, l)

    dera = derα*(8*sin(2α))/((cos(2α)-3)^2)
    derb = derα*(8*l*(cos(α)^3))/((cos(2α)-3)^2)
    derc = derb/(l^2)

    if wall=="horizontal"
        n1=0; n2=1;
    elseif wall=="vertical"
        n1=1; n2=0;
    end

    B_mat =
    [n2^2+a*n1^2 n1*n2*(a-1) -c*n1;
    n1*n2*(a-1) n1^2+a*n2^2 -c*n2;
    -b*n1 -b*n2 -a]

    A_mat =
    [0 0 n1*(derα*(px*n1+py*n2))-L*derc;
    0 0 n2*(derα*(px*n1+py*n2))-L*derc;
    0 0 -L*dera-derb*(n1*px+n2*py)]

    vecpp = A_mat*vecqq + B_mat*vecpp

    return [vecqq, vecpp]
end



function evol_function(Γ, l=l0, m=m0)
    I = m*l^2
    px, py, L = Γ[4:6]; vx = px/m; vy = py/m; ω=L/I

    [vx, vy, ω, 0, 0, 0]
end



function δτ(δΓ, Γ, part, wall, l=l0, m=m0)
    I = m*l^2
    dx, dy, dθ = δΓ[1:3]
    px, py, L = Γ[4:6]; vx = px/m; vy = py/m; ω = L/I
    θ=Γ[3]

    if part==1
        dxi = dx - l*sin(θ)*dθ
        dyi = dy + l*cos(θ)*dθ
        vxi = vx - l*ω*sin(θ)
        vyi = vy + l*ω*cos(θ)

    elseif part==2
        dxi = dx + l*sin(θ)*dθ
        dyi = dy - l*cos(θ)*dθ
        vxi = vx + l*ω*sin(θ)
        vyi = vy - l*ω*cos(θ)
    end

    dqqi = [dxi, dyi]
    vvi = [vxi, vyi]

    if wall=="vertical"
        n = [1,0]
    elseif wall=="horizontal"
        n = [0, 1]
    end

    return -1*(dqqi⋅n)/(vvi⋅n)
end


function displacement_vector_collision_map(Γ, δΓ, part, wall, l=l0, m=m0)
    Γp = collision_map(Γ, part, wall, l)
    DMdΓ = collision_map_derivative(Γ, δΓ, part, wall, l)
    FΓ = evol_function(Γ, l, m)
    DMF = collision_map_derivative(Γ, FΓ, part, wall, l)
    FΓp = evol_function(Γp, l ,m)

    dT = δτ(δΓ, Γ, part, wall, l, m)

    return DMdΓ + (DMF-FΓp)*dT
end


function evol_Γ(Γ0, t, l=l0, m=m0)
    I=m*l^2
    qq = Γ0[1:3]; pp = Γ0[4:6]
    vx, vy = pp[1:2]/m; ω = pp[end]/I
    qq += t*[vx, vy, ω]

    return [qq, pp]
end

function evol_δΓ(δΓ0, t, l=l0, m=m0)
    I = m*l^2
    dqq = δΓ0[1:3]; dpp = δΓ0[4:6]
    dvx, dvy = dpp[1:2]/m; dω = dpp[end]/I
    dqq += t*[dvx, dvy, dω]

    return [dqq, dpp]
end

function check_dumbbell_test(check_frequency, db::dumbbell, walls, exponents)
    cc = db.collision_counter
    cf = check_frequency
    if cc%cf==0
        exp_time, exp_cc = exponents
        sum1 = 0.; sum2 = 0.
        for j in 1:6
            sum1 += exp_time[cc,j]
            sum2 += exp_cc[cc,j]
        end

        println("Sum of Lyapunov Exponentes after \t"*string(cc)*" collisions")
        println(sum1, "\t", sum2)

        x1, y1 = particle_position(db, 1)[1:2];   x2, y2 = particle_position(db, 2)[1:2]
        vert_walls = walls[1:2]; horiz_walls = walls[3:4] # coordinates of the vertical and horizontal walls, respectively
        err = 1e-10

        if ~(minimum(vert_walls)-err< x1 < maximum(vert_walls)+err) ||
            ~(minimum(vert_walls)-err< x2 < maximum(vert_walls)+err) ||
            ~(minimum(horiz_walls)-err< y1 < maximum(horiz_walls)+err) ||
            ~(minimum(horiz_walls)-err< y2 < maximum(horiz_walls)+err)

            println("Dumbbell position\n(particle 1; particle 2)")
            println(x1,"\t", y1,";\t", x2,"\t", y2)

            error("The dumbbell is out of the billiard")

        end

    end

end


function disp_vecs_transformation(db::dumbbell, walls, norm_disp_vecs)

    ### Given a state of the dumbbell, 'db', inside thw billiard with frontiers, 'walls,
    ### And a set of normalized displacement vectors, 'norm_disp_vecs',
    ### This function will calculate the evolution of the displacement vectors (DV); their transformation
    ### during the collision, and a new set of orthogonal and orthonormal vectors: 'ortho_disp_vecs';
    ### norm_disp_vecs2.
    ### The output is a tuple of the collision time ('tc') it took the dumbbell to reach for the next wall,
    ### ortho_disp_vecs, and norm_disp_vecs2

    vert_walls = walls[1:2]; horiz_walls = walls[3:4] # coordinates of the vertical and horizontal walls, respectively

    Γ0,  Γf, wall, part, ct = collision(db, vert_walls, horiz_walls)

    disp_vecs0 = Array{Float64}[]
    for v in norm_disp_vecs
        vf = evol_δΓ(v, ct, l0, m0) # evolution of each of the initial normalized DV
        push!(disp_vecs0, [vf])
    end

    disp_vecs_F = Array{Float64}[]
    for dΓ in disp_vecs0
        dΓp = displacement_vector_collision_map(Γ0, dΓ, part, wall, l0, m0)
        push!(disp_vecs_F, [dΓp])
    end

    ortho_disp_vecs, norm_disp_vecs2 = Gram_Schmidt(disp_vecs_F)

    ct, ortho_disp_vecs, norm_disp_vecs2
end

function lyapunov_spectrum(db::dumbbell, walls, N, check=false, check_times=3)
    # N=number of collisions
    srand(5) #seed for random number generator

    # in case 'check' variable is True, every 'check_frequency' collisions,
    # the sum of the exponents will be computed and a test will be run to guarantee that the dumbbell
    # is within the frontiers of the billiard (+/- the error with which a collision is determined)
    check_frequency = int(N/check_times)

    rand_vecs = ([rand(6) for i in 1:6]) # initial random displacement vectors
    ortho_disp_vecs, norm_disp_vecs = Gram_Schmidt(rand_vecs)  # initial DV, orthogonalized and orthonormalized

    col_count = 0 # collision counter

    λ_col = zeros(Float64,N,6); λ_time = zeros(Float64,N,6) #value of the LEs
    χs = zeros(Float64,N,6) # log of the norm of (the sum of) the DV

    # arrays of the collision time and the total elapsed time (until each collision)
    col_times = Float64[]; elapsed_time = Float64[];
    t_T= 0. # initial time

    # main loop for calculating the LEs for the N collisions
    for k in 1:N
        col_count += 1

        #collision time, orthogonal vectors, orthonormal vectors after the k-th collision
        ct, O_vecs, ON_vecs = disp_vecs_transformation(db, walls, norm_disp_vecs)

        t_T += ct # add the current collision time to the elapsed time
        push!(elapsed_time, t_T);  push!(col_times, ct) # register the times

        for (j,ov) in enumerate(O_vecs)
            if k>1
                χs[k, j] = χs[k-1, j] + log(norm(ov))
            else
                χs[1,j] = log(norm(ov))
            end

            # LEs obtained with the number of collisions and the elapsed time
            λ_col[k,j] = χs[k,j]/(k); λ_time[k,j] = χs[k,j]/t_T
        end

        if check
            check_dumbbell_test(check_frequency, db, walls, (λ_col, λ_time))
        end


        norm_disp_vecs = ON_vecs # set the new DV as the normalized DV after the collision

    end

    λ_col, λ_time, elapsed_time, col_times
end

