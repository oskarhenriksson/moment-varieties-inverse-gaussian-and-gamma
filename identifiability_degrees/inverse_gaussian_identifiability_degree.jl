using HomotopyContinuation
using Random
include("certify_orbits.jl")

# Set seed for reproducibility
Random.seed!(12345)

k = 2;
d = 3*k-1;

# Define the moments polynomials
@var mu, kappa
moment_polynomials = [mu,expand(mu^2*(1+mu*kappa))]
for i = 3:d
    append!(moment_polynomials, expand((2*i-3)*kappa*mu^2*moment_polynomials[i-1] + mu^2*moment_polynomials[i-2]))
end

# Define the method of moments system
@var α[1:k] μ[1:k] κ[1:k] m[1:d] a b[1:k]
parametrization = vcat([sum([α[i] for i=1:k])],
        [sum([α[i]*subs(moment_polynomials[r],[mu=>μ[i],kappa=>κ[i]]...) for i=1:k]) for r=1:d ])
F = System(parametrization-vcat(a,m),
        parameters=vcat(a,m),
        variables=vcat(α,κ,μ)
    )

# Define the label-swapping action
relabeling = v -> map(p -> [v[1:k][p]..., v[k+1:2*k][p]..., v[2*k+1:3*k][p]...], SymmetricGroup(k))

# Solve the system for random sample moments
random_distribution_parameters = rand(Complex{Float64},3*k)
random_system_parameters = parametrization(F.variables=>random_distribution_parameters)
R = monodromy_solve(F,[random_distribution_parameters],random_system_parameters, group_action=relabeling, seed=0x12345)

# Certify the result
certify_orbits(F, R, relabeling)