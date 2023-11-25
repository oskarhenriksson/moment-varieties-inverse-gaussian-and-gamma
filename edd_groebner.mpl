# This file contains code for the proof of Proposition 4.5


# Degree reverse lexicographic order
ord := tdeg(seq(m[i],i=1..d),a,b)

# Gaussian := -2*m[1]^3 + 3*m[1]*m[2] - m[3];
G := Groebner[Basis]([f,seq( diff(f,m[i])-a*(m[i]-n[i]), i=1..3 ),a*b-1],ord):
Groebner[NormalSet](G,ord)[1]; 
nops(%);

# Gamma := -m[1]*m[3] + 2*m[2]^2 - m[1]^2*m[2];
G := Groebner[Basis]([f,seq( diff(f,m[i])-a*(m[i]-n[i]), i=1..3 ),a*b-1],ord):
Groebner[NormalSet](G,ord)[1]; 
nops(%);

# Inverse Gaussian
f := -m[1]*m[3] + 3*m[2]^2 - 3*m[1]^2*m[2] + m[1]^4;
G := Groebner[Basis]([f,seq( diff(f,m[i])-a*(m[i]-n[i]), i=1..3 ),a*b-1],ord):
Groebner[NormalSet](G,ord)[1]; 
nops(%);
