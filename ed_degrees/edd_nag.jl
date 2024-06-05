using HomotopyContinuation

@var m[1:d] mhat[1:d] a b

list_of_polynomials = [-2*m[1]^3 + 3*m[1]*m[2] - m[3], #gaussian
    -m[1]*m[3] + 2*m[2]^2 - m[1]^2*m[2], #gamma
    -m[1]*m[3] + 3*m[2]^2 - 3*m[1]^2*m[2] + m[1]^4, #inverse gaussian
    ];

for f in list_of_polynomials
    
    println("Polynomial defining hypersurface: $(f)\n")
    
    F = System(vcat([f], [differentiate(f,m[i])-a*(m[i]-mhat[i]) for i=1:d], [a*b-1]),
    variables=vcat(m,[a,b]),parameters=mhat)
    
    println("Mixed volume: $(mixed_volume(F))\n")
    
    sample_moments = rand(Complex{Float64},d)
    S = solve(F,target_parameters=sample_moments)
    display(S)
    println()
    
    C = certify(F,S,target_parameters=sample_moments)
    display(C)
    println()
    
end
