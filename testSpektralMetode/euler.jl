using Plots, LinearAlgebra

function euler_matrix(N::Int, α::Float64, β::Float64, r::Vector{Float64})


    beta_array = β*ones(N)./r


    diag = ones(N-1)*(1 -2* α)

    off_diagup = ones(N-2)*(α)
    off_diaglo = ones(N-2)*(α)



   off_diaglo .-= beta_array[3:end]

   off_diagup .+= beta_array[2:end-1]

    A = Tridiagonal(off_diaglo, diag, off_diagup)


    return A


end


function initialcondition(r::Vector{Float64})


    return  r.^2 .+4

end

function initialcond2(N::Int)

    u1 = zeros(N)

    for i in 1:N
        u1[i] = 5 + (i-5)*(i-100)/1000

    return u1
    end

end




origoInterpolasjon(u) = (1/3)*(4*u[2] - u[3])


function initialilize(N::Int, J::Int, d_t::Float64, D::Float64, T::Int)


    r = range(0, 1, length=N)

    Θ = range(0, 2*π, length=J)

    d_r = r[2] - r[1]

    α = (D*d_t)/(d_r^2)

    β = (D*d_t)/(2*d_r)


    A = euler_matrix(N, α, β, r)


    A[1,2] = 2*α

    A[end, end] = 1

    A[end, end-1] = 0

    
    u = zeros(N, T)

    return r, Θ, A, u

end


function run(N::Int, J::Int, T::Int, d_t::Float64, D::Float64)


    r, Θ, A, u = initialilize(N, J, d_t, D, T)

    u[:, 1] = initialcondition(r)


    for i in 1:T-1

        u[2:end, i+1] = A*u[2:end, i]
        

        u[1, i+1] = origoInterpolasjon(u[:, i+1])


    end


    return u, A, r

end

u, A, r = run(100, 10, 100000, 1E-5, 1)


plot(u[:, 1])
plot!(u[:, Int(end/2)])