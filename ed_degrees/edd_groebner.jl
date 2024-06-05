# This file contains code for the proof of Proposition 4.5

using Oscar

# Set up the coefficient field and ambient ring
d = 3
K, m̂ = rational_function_field(QQ, "mhat"=>1:d)
R, m, k = polynomial_ring(K, "m"=>1:d, "k"=>1:1)
k = k[1]

# Inverse Gaussian
f = -m[1]*m[3] + 3*m[2]^2 - 3*m[1]^2*m[2] + m[1]^4
I = ideal(vcat([f],[k*derivative(f,m[i])-(m[i]-m̂[i]) for i=1:3]))
R_mod_I, _ = quo(R,I)
println(vector_space_dimension(R_mod_I))

# Gamma 
f = -m[1]*m[3] + 2*m[2]^2 - m[1]^2*m[2];
I = ideal(vcat([f],[k*derivative(f,m[i])-(m[i]-m̂[i]) for i=1:3]))
R_mod_I, _ = quo(R,I)
println(vector_space_dimension(R_mod_I))

# Gaussian 
f = -2*m[1]^3 + 3*m[1]*m[2] - m[3];
I = ideal(vcat([f],[k*derivative(f,m[i])-(m[i]-m̂[i]) for i=1:3]))
R_mod_I, _ = quo(R,I)
println(vector_space_dimension(R_mod_I))
