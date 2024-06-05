using HomotopyContinuation

k = 3;
d = 3*k-1;

@var α[1:k] μ[1:k] λ[1:k] κ[1:k] m[1:d] a b[1:k]

@var mu, lambda, kappa

moment_polynomials = [mu,mu^2*(lambda+mu)*kappa]
for i = 3:d
	append!(moment_polynomials, expand((2*i-3)*kappa*mu^2*moment_polynomials[i-1] + mu^2*moment_polynomials[i-2]))
end

parametrization = vcat([sum([α[i] for i=1:k])],
        [λ[i]*κ[i] for i=1:k],
        [sum([α[i]*subs(moment_polynomials[r],[mu=>μ[i],lambda=>λ[i],kappa=>κ[i]]...) for i=1:k]) for r=1:d ])

F = System(parametrization-vcat(a,b,m),
        parameters=vcat(a,b,m),
        variables=vcat(α,κ,λ,μ)
    )

G = SymmetricGroup(k)
relabeling(v) = map(p -> (v[1:k][p]..., v[k+1:2*k][p]..., v[2*k+1:3*k][p]..., v[3*k+1:4*k][p]...), G)


random_distribution_parameters = rand(Complex{Float64},4*k)

random_system_parameters = parametrization(F.variables=>random_distribution_parameters)

R = monodromy_solve(F,[random_distribution_parameters],random_system_parameters, group_action=relabeling, seed=0x42953)

C = certify(F,R)
