using Integrals, LinearAlgebra, CairoMakie, FFTW


function fourierseries(N, T, Lx, nx)

    x1 = LinRange(-Lx/2, Lx/2, nx)

    y2 = zeros(nx)


    for i in 1:N

        y2 += (cos.(((2 * π * i)/Lx)*x1) + sin.(((2 * π * i)/Lx)*x1 ))

    end


return y2, x1

end


function integraltesst()

    x = LinRange(-0.5, 0.5, 1000000)

    y = 2*(0.2*cos.(2*10*π*x))

    metod = TrapezoidalRule()

    problem = SampledIntegralProblem(y, x)

    return solve(problem, metod)
end



function fftwtest()

    y = LinRange(-0.5, 0.5, 20)



    ff1 =fft(0.2*sin.(((2 * π )*y )))

    a = real(ff1)[1:10]

    println(a)

    b = imag(ff1)

    ux = zeros(20)

    for i in 1:10
        
        ux += a[i]*cos.(((2 * π * i)/1)*y) + b[i]*sin.(((2 * π * i)/1)*y) 


    end

    lines(y, (1/10)*ux)
end

fftwtest()