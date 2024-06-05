using HomotopyContinuation

k = 3;
d = 3*k-1;

@var α[1:k] β[1:k] λ[1:k] m[1:d];

parametrization = vcat([sum([λ[i] for i=1:k])],
	[ sum([λ[i]*β[i]^r * prod(Expression[(α[i]+j) for j=1:(r-1)]) for i=1:k]) for r=1:d ]);

F = System(parametrization-vcat([1],m),parameters=m)

G = SymmetricGroup(k);
relabeling(v) = map(p -> (v[1:k][p]..., v[k+1:2*k][p]..., v[2*k+1:3*k][p]...), G);

R = monodromy_solve(F,group_action=relabeling,seed=0x42953765)

C = certify(F,R)
