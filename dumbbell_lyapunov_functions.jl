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