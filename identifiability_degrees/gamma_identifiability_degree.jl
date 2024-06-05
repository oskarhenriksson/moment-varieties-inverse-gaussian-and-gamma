using HomotopyContinuation
using Random
include("certify_orbits.jl")

# Set seed for reproducibility
Random.seed!(123456)

k = 2
d = 3*k-1;

# Define the method of moments system
@var α[1:k] β[1:k] λ[1:k] m[1:d];
parametrization = vcat([sum([λ[i] for i=1:k])],
	[ sum([λ[i]*β[i]^r * prod(Expression[(α[i]+j) for j=0:(r-1)]) for i=1:k]) for r=1:d ]);
F = System(parametrization-vcat([1],m),parameters=m)

# Define the label swapping symmetry
G = SymmetricGroup(k);
relabeling(v) = map(p -> [v[1:k][p]..., v[k+1:2*k][p]..., v[2*k+1:3*k][p]...], G);

# Solve the system
R = monodromy_solve(F,group_action=relabeling, max_loops_no_progress=10, seed=0x123456)

# Certify the result
certify_orbits(F, R, relabeling, save_results=true)