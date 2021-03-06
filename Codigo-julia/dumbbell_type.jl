import Base.show
type dumbbell
    position::Vector
    angle::Float64
    velocity::Vector
    omega::Float64
    m::Float64
    l::Float64
    collision_counter::Int64
    
end

show(io::IO, db::dumbbell)= print(io, "CM: $(db.position); θ: $(round((db.angle)%(2*pi),2)); "*
"Vcm: $(db.velocity); ω: $(db.omega); m=$(db.m); l=$(db.l); I: $(round(db.m*db.l^2,2))")


function move(db::dumbbell, dt)
    db.position += db.velocity*dt
    db.angle += db.omega*dt
    #db.angle = db.angle%(2pi)
end


function particle_position(db::dumbbell, part::Int, t=0) # define la posición de la partícula 'part' al tiempo 't'
    x, y = db.position
    vx, vy = db.velocity
    θ = db.angle;    ω = db.omega
    l = db.l
    θ_f = (θ + ω*t)%(2pi) 
    
    if part==1
        x1 = x + vx*t + l*cos(θ_f)
        y1 = y + vy*t + l*sin(θ_f) 
        
        vx1 = vx - l*ω*sin(θ_f)
        vy1 = vy + l*ω*cos(θ_f)
        return [x1, y1, vx1, vy1]
        
    elseif part==2
        x2 = x + vx*t - l*cos(θ_f)
        y2 = y + vy*t - l*sin(θ_f) 
        
        vx2 = vx + l*ω*sin(θ_f)
        vy2 = vy - l*ω*cos(θ_f)
        
        return [x2, y2, vx2, vy2]
        
    else
        println("No hay tantas partículas")
    end
end


function reset_dummbell(db::dumbbell, original_position = [1., 1.], 
    original_θ=pi/6, original_velocity = [2.1, 0.6],
    original_ω=2.0)
    db.position = original_position
    db.velocity = original_velocity
    db.angle = original_θ
    db.omega = original_ω
    db.collision_counter = 0;
    
    return "Dumbbell reset"
end


function ct_approx(db::dumbbell, wall, direction)
    # This function returns an *approximation* of the time needed for one of the particles of the dumbbell
    # to collide with the walls located at 'w1' and 'w2'. 'direction' specifies the orientation of the walls
    #involved: if they are vertical, only the motion on the x axis is relevant, hence, direction=1;
    #when they are horizontal, only the motion on the y axis needs to be considered, therefore, direction=2 
    #dt=0.01     
    err = 1e-10
    t_elapsed = 0.
    err2 = 0.01
    w1 = wall[1]; w2 = wall[2]
    
    v = db.velocity; l = db.l; ω=db.omega
    
    # positions of the particles 1 and 2, respectively, along the specified direction
    r1 = particle_position(db, 1)[direction]
    r2 = particle_position(db, 2)[direction]
    
    # distances from the both particles to each wall    
    d11 = abs(w1-r1); d12 = abs(w1-r2); d21 = abs(w2-r1); d22 =abs(w2-r2)
    V = sqrt(norm(v)^2+ (l*ω)^2)
    
    if db.collision_counter!=0
        distances = sort([d11, d12, d21, d22])
        dmin = distances[2]
        times = [dmin, 1.]/V
        
    else
        distances = [d11, d12, d21, d22, 1.]
        times = distances/V
    end

    
    dt = 0.1*minimum(times)
    
    while w1-err <= r1 <= w2+err && w1-err <= r2 <= w2+err
        
        t_elapsed += dt
            
        r1 = particle_position(db, 1, t_elapsed)[direction]
        r2 = particle_position(db, 2, t_elapsed)[direction]
        
    end
    # The output is the last time when the condition of both particles being within 'w1' y 'w2'  was fulfilled;
    #and the time step used for the approximation  
    t_elapsed - dt, dt 
    
end

# this function gives the inertia moment of a dummbell with TOTAL mass 'm' and length '2l'
inertia_moment(m=1., l=0.2) = m*l^2  

function collision_transformation(v_i, ω_i, α, l=0.2, m=1., I = inertia_moment(1, 0.2))
    # with this function, the transformation of the velocity (normal to the wall)
    # and the angular velocity after a collision is obtained; 
    # the input needed are the initial velocities (v_i, ω_i),
    # the angle between the angular coordinate of the dumbbell (θ) and the normal to the wall of collission, α;
    # and the mass, length and intertia moment of the dumbbell.
    # The output are the new (normal component of the) velocity and angular velocity, v_f, ω_f
    
    denom = m*l^2 + I*csc(α)^2    
       
    v_f = (1 - (2*I*csc(α)^2)/denom)*v_i - (2*l*I*csc(α)/denom)*ω_i
    
    ω_f = (-2*l*m*csc(α)/denom)*v_i - (1-(2*I*csc(α)^2)/denom)*ω_i
    
    
    v_f, ω_f
end


function col_condition(db::dumbbell, i::Int, j::Int, t, wall, direction::Int, t_approx, δt=0.01)
    err = 1e-10
    w = wall[j]
    r = particle_position(db, i, t)[direction]
        
    dist = abs(w-r);
   (dist<err && t_approx < t <= t_approx + δt)
end


function directional_collision(db::dumbbell, walls, axis)
    vx, vy = db.velocity
    θ, ω, l, m = db.angle, db.omega, db.l, db.m
    I = inertia_moment(m, l)
    err = 1e-10
               
    tapprox, dt = ct_approx(db, walls, axis)
    times = Float64[]
    max_iter = 20
    
    col_times = zeros(Float64, 2, 2)
    position = zeros(Float64, 2, 2)
    velocity = zeros(Float64, 2, 2)
    counter = zeros(Int32, 2 ,2)
    
    for i in [1:2], j in [1:2]
        col_times[i,j] = tapprox
        position[i,j], velocity[i,j] = particle_position(db, i, col_times[i,j])[[axis, axis+2]]
    end
    
    for i in [1:2], j in [1:2]
        while ~col_condition(db, i, j, col_times[i,j], walls, axis, tapprox, dt)
            col_times[i,j] = col_times[i,j] - (position[i,j]-walls[j])/velocity[i,j]
            position[i,j], velocity[i,j] = particle_position(db, i, col_times[i,j])[[axis, axis+2]]
            counter[i,j] += 1
            if counter[i,j]>max_iter
                col_times[i,j] = -1.
                break
            end
        end
        push!(times, col_times[i,j])
    end
         
   # elegir cuál es el tiempo correcto
    times = times[times.>0]
    if length(times)!=0
        tc = minimum(times)
    else
        return("Error  ", db.collision_counter, col_times)
    end
           
    θ = (θ + ω*tc)%(2pi) # angle position of the dumbbell when the collision occurs

    if axis ==1  # collision in the 'x' direction (vertical walls)
        if tc==col_times[1, 1] || tc==col_times[1, 2] # the first particle collides
            α = -θ
            part = 1
        elseif tc==col_times[2, 1] || tc==col_times[2, 2] # the second particle collides
            α = θ
            part = 2
        end
        
        v_f, ω_f = collision_transformation(vx, ω, α, l, m, I)
        
        return tc, v_f, vy, ω_f, part, axis, col_times
        
        
    elseif axis==2 # collision in the y direction (horizonta walls)
        if tc==col_times[1, 1] || tc==col_times[1, 2] # the first particle collides
            α = pi/2 - θ
            part = 1
        elseif tc==col_times[2, 1] || tc==col_times[2, 2] # the second particle collides
            α = θ - pi/2
            part = 2
        end
        
        v_f, ω_f = collision_transformation(vy, ω, α, l, m, I)
        
        return tc, vx, v_f, ω_f, part, axis, col_times
        
    else
        println("Incorrect direction")
    end
            
end



function collision(db::dumbbell, vert_walls = [-3, 3], horiz_walls =[-2, 2])
    db.collision_counter += 1
    tcx, vx_fx, vy_fx, ω_fx, part_x = directional_collision(db, vert_walls, 1)[1:5]
    tcy, vx_fy, vy_fy, ω_fy, part_y = directional_collision(db, horiz_walls, 2) [1:5]
    
    x_0, y_0 = db.position; θ_0 = db.angle
    vx_0, vy_0 = db.velocity; ω_0 = db.omega
    l = db.l; m = db.m; I = inertia_moment(m,l)
    
    # vector del espacio fase inicial (antes de que la mancuerna se mueva hacia el punto de colisión)
    Γ_0 = [x_0, y_0, θ_0, vx_0, vy_0, I*ω_0]
    
    if tcx < tcy && tcx >0
        tc, vx_new, vy_new, ω_new, part = tcx, vx_fx, vy_fx, ω_fx, part_x
        wall = "vertical"
        
    elseif tcy < tcx && tcy>0
        tc, vx_new, vy_new, ω_new, part = tcy, vx_fy, vy_fy, ω_fy, part_y
        wall = "horizontal"
    else
        return "Error"
        
    end
    
    move(db,tc) # the dumbbell is moved to the collision point
           
    db.velocity = [vx_new, vy_new]
    db.omega = ω_new
    x_f, y_f = db.position; θ_f = db.angle
    
    Γ_f = [x_f, y_f, θ_f, vx_new, vy_new, I*ω_new]
    
    return (Γ_0, Γ_f, wall, part, tc)
    
end

function next_collision(db::dumbbell, vert_walls = [-3, 3], horiz_walls =[-2, 2])
    
    tcx, vx_fx, vy_fx, ω_fx, part_x = directional_collision(db, vert_walls, 1)[1:5]
    tcy, vx_fy, vy_fy, ω_fy, part_y = directional_collision(db, horiz_walls, 2)[1:5]
    
    x_0, y_0 = db.position; θ_0 = db.angle
    vx_0, vy_0 = db.velocity; ω_0 = db.omega
    l = db.l; m = db.m; I = inertia_moment(m,l)
    
    if tcx < tcy && tcx >0
        tc, vx_new, vy_new, ω_new, part = tcx, vx_fx, vy_fx, ω_fx, part_x
        wall = "vertical"
        
    elseif tcy < tcx && tcy>0
        tc, vx_new, vy_new, ω_new, part = tcy, vx_fy, vy_fy, ω_fy, part_y
        wall = "horizontal"
    else
        return "Error"
    end
        
    x1, y1 = particle_position(db, 1, tc)
    x2, y2 = particle_position(db, 2, tc)
    
    new_r = db.position + db.velocity*tc
    new_θ = (db.angle + db.omega*tc)%(2pi)
    new_v = [vx_new, vy_new]
    new_ω = ω_new

    new_db = dumbbell(new_r, new_θ, new_v, new_ω, db.m, db.l, 0)
    
    return new_db, [tc, x1, y1, x2, y2]

end


function collision_time(db::dumbbell, vert_walls=[-3,3], horiz_walls=[-2,2])
    tcx = directional_collision(db, vert_walls, 1)[1]
    tcy = directional_collision(db, horiz_walls, 2)[1]
    if tcx<tcy
        return [tcx, 1]
    else
        return [tcy, 2]
    end
end



function collision_parts(db::dumbbell, vert_walls=[-3,3], horiz_walls=[-2,2])
    db.collision_counter +=1
    ct, axis = collision_time(db, vert_walls, horiz_walls)
    
    tcx, vx_fx, vy_fx, ω_fx, part_x = directional_collision(db, vert_walls, 1)[1:5]
    tcy, vx_fy, vy_fy, ω_fy, part_y = directional_collision(db, horiz_walls, 2) [1:5]
    
    x_0, y_0 = db.position; θ_0 = db.angle
    vx_0, vy_0 = db.velocity; ω_0 = db.omega
    l = db.l; m = db.m; I = inertia_moment(m,l)
  
    
    if tcx < tcy && tcx >0
        tc, vx_new, vy_new, ω_new, part = tcx, vx_fx, vy_fx, ω_fx, part_x
        wall = "vertical"
        
    elseif tcy < tcx && tcy>0
        tc, vx_new, vy_new, ω_new, part = tcy, vx_fy, vy_fy, ω_fy, part_y
        wall = "horizontal"
    else
        return "Error"
        
    end
    
    move(db,tc) # the dumbbell is moved to the collision point
           
    db.velocity = [vx_new, vy_new]
    db.omega = ω_new
    
    x1, y1 = particle_position(db,1)[1:2]
    x2, y2 = particle_position(db,2)[1:2]
    
    return [x1, y1, x2, y2, ct, wall]
    
end
    

function orbit(db::dumbbell, T, dt=0.05, vert_walls=[-3,3], horiz_walls=[-2,2])

    ct = collision_time(db, vert_walls, horiz_walls)[1]
    t = 0.
    x1, y1 = particle_position(db,1, t)[1:2]
    x2, y2 = particle_position(db,2,t)[1:2]
    pos_p1 = [x1 y1]
    pos_p2 = [x2 y2]    
    times = [t]
    
    # T is less than the collisiont time 'ct'
    if T<ct
        n = int(T/dt)
        for t in linspace(dt, T, n)
            x1, y1 = particle_position(db,1, t)[1:2]
            x2, y2 = particle_position(db,2,t)[1:2]
            pos_p1 = vcat(pos_p1, [x1 y1])
            pos_p2 = vcat(pos_p2, [x2 y2])
            push!(times, t)
        end
        
        return times, pos_p1, pos_p2
    
        # T is greater than the collision time
    else
       while t <= T
            nt=0.
            while nt < ct && t<T
                x1, y1 = particle_position(db,1, nt)[1:2]
                x2, y2 = particle_position(db,2,nt)[1:2]
                pos_p1 = vcat(pos_p1, [x1 y1])
                pos_p2 = vcat(pos_p2, [x2 y2])
                push!(times, t)
                t += dt
                nt += dt
            end
            collision_parts(db, vert_walls, horiz_walls)
            ct = collision_time(db, vert_walls, horiz_walls)[1]
        end
        return times, pos_p1, pos_p2
    end
    
end

