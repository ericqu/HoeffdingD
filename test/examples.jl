using Plots# , HoeffdingD
using StatsBase
using BenchmarkTools

# ρ(N, α) = √(2 * (N^2 + 5 * N - 32) / (9N * (N - 1) * (N - 3) * (N - 4) * α ))

# n = 5:100
# αs = 0:0.02:1
# ys = ρ.(n, α)

# plot(n, ρ.(n, 0.05), label="α=0.05", xlabel="n")
# plot!(n, ρ.(n, 0.1), label="α=0.1")
# plot!(n, ρ.(n, 0.2), label="α=0.2")
# plot!(n, ρ.(n, 0.02), label="α=0.5")
# plot!(n, ρ.(n, 0.01), label="α=0.5")


# @btime corspearman(x, y_linear)
# @btime StatsBase.cor(x, y_linear)
# @btime HoeffdingD.hoeffdingd(x, y_linear)

# using Distributions, Distances, EnergyStatistics
using Distributions, HoeffdingD

function compareD(interval=-√(0.99):0.05:√(0.99))
    ρ = collect(interval) # -1:0.02:1
    meanPearson = 0.
    meanHoeffding = 0.

    Pmeans = Vector{Float64}(undef, length(ρ))
    Hmeans = Vector{Float64}(undef, length(ρ))
    nb_sims = 10
    for i in 1:length(ρ)
        meanPearson = 0.
        meanHoeffding = 0.
        # meanDcor = 0.
        meansP = Vector{Float64}(undef, nb_sims)
        meansH = Vector{Float64}(undef, nb_sims)
        for j in 1:nb_sims
            z = rand(MvNormal([0., 0.], ρ[i]), 500)
            a = view(z, 1, :)
            b = view(z, 2, :)
            meansP[j] = StatsBase.cor(a, b)
            meansH[j] = first(HoeffdingD.hoeffdingd(a, b))
        end
        Pmeans[i] = mean(meansP)
        Hmeans[i] = mean(meansH)
    end

    return Pmeans, Hmeans, ρ
end



plot(ρ, Pmeans)
plot!(ρ, Hmeans)

# N = 1000;
# Mean = {1 2};
# Cov = {2.4 3, 3 8.1};

# x = RandNormal( N, Mean, Cov );
# SampleMean = mean(x);
# SampleCov = cov(x);
# print SampleMean Mean, SampleCov Cov;

covsigma = 
@show(covsigma)

MvNormal([0, 0], [-1, -1])

# z1 = rand(MvNormal([0, 0], [-1 0 ; 0 -1]), 1000)
z1 = rand(MvNormal([2, 1], [-1, -1]), 20000)
# z1 = rand(MvNormal([2, 1], [2.4 1; 1 8.3 ]), 20000)
a = view(z1, 1, :)
b = view(z1, 2, :)
@show(mean(a))
@show(mean(b))
@show(cov(a))
@show(cov(b))
@show(StatsBase.cor(a, b))
@show(HoeffdingD.hoeffdingd(a, b))


# The CORR Procedure

# 6 Variables:	MSRP Invoice EngineSize Horsepower MPG_City MPG_Highway
# Pearson Correlation Coefficients, N = 428
#  	MSRP	Invoice	EngineSize	Horsepower	MPG_City	MPG_Highway
# MSRP	1.00000	0.99913	0.57175	0.82695	-0.47502	-0.43962
# Invoice	0.99913	1.00000	0.56450	0.82375	-0.47044	-0.43459
# EngineSize	0.57175	0.56450	1.00000	0.78743	-0.70947	-0.71730
# Horsepower	0.82695	0.82375	0.78743	1.00000	-0.67670	-0.64720
# MPG_City	-0.47502	-0.47044	-0.70947	-0.67670	1.00000	0.94102
# MPG_Highway	-0.43962	-0.43459	-0.71730	-0.64720	0.94102	1.00000
# Hoeffding Dependence Coefficients, N = 428
#  	MSRP	Invoice	EngineSize	Horsepower	MPG_City	MPG_Highway
# MSRP	0.99945	0.93253	0.16067	0.35767	0.19844	0.14221
# Invoice	0.93253	0.99995	0.15670	0.35327	0.19646	0.14201
# EngineSize	0.16067	0.15670	0.88251	0.31255	0.35381	0.24217
# Horsepower	0.35767	0.35327	0.31255	0.95970	0.27758	0.19187
# MPG_City	0.19844	0.19646	0.35381	0.27758	0.76946	0.50466
# MPG_Highway	0.14221	0.14201	0.24217	0.19187	0.50466	0.82463


