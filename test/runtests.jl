using Test, HoeffdingD
ab =[
    -1.0	2.03
    -1.0	2.11
    -1.0	2.03
    -0.8	1.58
    -0.8	1.79
    -0.8	1.58
    -0.6	1.35
    -0.6	1.46
    -0.6	1.37
    -0.4	1.28
    -0.4	1.15
    -0.4	1.19
    -0.2	0.93
    -0.2	0.99
    -0.2	0.99
    0.0 1.03
    0.0 0.98
    0.0 1.02
    0.2 1.00
    0.2 1.06
    0.2 1.00
    0.4 1.31
    0.4 1.33
    0.4 0.99
    0.6 1.28
    0.6 1.42
    0.6 1.50
    0.8 1.68
    0.8 1.64
    0.8 1.68
    1.0 2.07
    1.0 1.94
    1.0 2.14 ]

a = ab[:, 1]
b = ab[:, 2]


@testset "HoeffdingD" begin
    #results values are from SAS PROC Corr 
    @test isapprox(0.1108736664, HoeffdingD.hoeffdingd(a, b))
    @test isapprox(0.1108736664, HoeffdingD.hoeffdingd(b, a))
    @test isapprox(0.9665923259, HoeffdingD.hoeffdingd(b, b))
    @test isapprox(0.8450535317, HoeffdingD.hoeffdingd(a, a))

    @test isapprox(0.6875, HoeffdingD.hoeffdingd(
            [9, 8,  7,  missing,    missing,    2,          5, 10, 2], 
            [6, 5,  4,  missing,    3,          missing,    -5, 1000, 4]))

    @test isapprox([0.845053531701891 0.110873666447568; 0.110873666447568 0.9665923258586983], HoeffdingD.hoeffdingd(ab))


    linear_f(x) = @. 3x + 2
    quad_f(x) = @. x^2 + 1
    x_test = -1:0.05:1
    y_linear = linear_f.(x_test)
    y_quad = quad_f.(x_test)
    xyy = hcat(x_test, y_linear, y_quad)
    @test isapprox( [1.0 1.0 0.2183712793468891; 1.0 1.0 0.2183712793468891; 0.2183712793468891 0.2183712793468891 0.9369001852153328] , HoeffdingD.hoeffdingd(xyy))
        
end

@testset "HoeffdingD α" begin
    @test_throws DomainError HoeffdingD.hoeffdingd(ab, -1)
    @test [true false; false true] == last.(HoeffdingD.hoeffdingd(ab, .001))
    @test [true true; true true] == last.(HoeffdingD.hoeffdingd(ab, .05))
end

using Distributions

@testset "HoeffdingD independence" begin
    x = collect(1:200)
    y_uniform = rand(Distributions.Uniform(2,5), length(x))
    @test isapprox(0, last(HoeffdingD.hoeffdingd(x, y_uniform, 0.95)))
    @test isapprox(0, last(HoeffdingD.hoeffdingd(x, y_uniform, 0.5)))
    @test isapprox(0, last(HoeffdingD.hoeffdingd(x, y_uniform, 0.1)))
    @test isapprox(0, last(HoeffdingD.hoeffdingd(x, y_uniform, 0.05)))
    @test isapprox(0, last(HoeffdingD.hoeffdingd(x, y_uniform, 0.005)))
    @test isapprox(0, last(HoeffdingD.hoeffdingd(x, y_uniform, 0.001)))
end


using BenchmarkTools

@btime HoeffdingD.hoeffdingd(a, b)
@btime HoeffdingD.hoeffdingd(b, a)
@btime HoeffdingD.hoeffdingd(ab)

@btime HoeffdingD.hoeffdingd(a, b, 0.05)
@btime HoeffdingD.hoeffdingd(b, a, 0.05)
@btime HoeffdingD.hoeffdingd(ab, 0.05)

@btime HoeffdingD.hoeffdingd( [9, 8,  7.,  missing,    missing,    2,          5, 10, 2], 
            [6, 5,  4.,  missing,    3,          missing,    -5, 1000, 4])

