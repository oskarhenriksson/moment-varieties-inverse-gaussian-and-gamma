restart
needsPackage "Depth"

--dth moment variety of the inverse Gaussian
d = 5
S = QQ[m_0..m_d, x, y, z]

moments = {1,y}

for i from 2 to d do (
    f = (2*i-3)*z*y^2*m_(i-1) + y^2 * m_(i-2);
    moments = append(moments,f);
    )

genList = for i from 0 to d list (
    m_i - moments_i
    )
genList = append(genList, x*z-1)

I = ideal genList
Ie = eliminate({x,y,z},I)
Ieh = homogenize(Ie, m_0)
Iehs = saturate(Ieh, m_0)


R = QQ[m_0..m_d]

--J is the defining ideal of M_d
J= sub(Iehs,R)
netList first entries gens J

dim J
degree J
hilbertSeries(J, Reduce=>true)
inJ= ideal mingens ideal leadTerm gb J

R1 = R/J;
isCM R1

