# Fix seed for reproducibility
randomize(12345);

for d from 1 to 100 do
	for k from 1 to d do

		print((d,k));

		# Expected dimension
		expected_dimension := min([d,2*k + k - 1]);

		moment_polynomials := [mu,mu^2*(lambda+mu)/lambda]:
		for i from 3 to d do
			moment_polynomials := [op(moment_polynomials), expand((2*i-3)/lambda*mu^2*moment_polynomials[i-1] + mu^2*moment_polynomials[i-2]) ]:
		end do:

		# Form the generators of the ideal of the incidence variety
		parametrization := subs([alpha[k]=1-add(alpha[i],i=1..(k-1))],[
			seq(add(alpha[i]*subs([mu=mu[i],lambda=lambda[i]],moment_polynomials[r]),i=1..k),r=1..d)]):

		Jac := VectorCalculus[Jacobian](parametrization,[seq(lambda[i],i=1..k),seq(mu[i],i=1..k),seq(alpha[i],i=1..(k-1))]):

		# Compute the rank at a random integer point
		subs([seq(lambda[i]=rand(1..1000)(),i=1..k),seq(mu[i]=rand(-1000..1000)(),i=1..k),seq(alpha[i]=rand(-1000..1000)(),i=1..(k-1))],Jac):
		lower_bound_of_dimension := LinearAlgebra[Rank](%);

		# Pick another point if the rank is lower than the expected one
		if lower_bound_of_dimension <> expected_dimension then
			subs([seq(lambda[i]=rand(1..10000)(),i=1..k),seq(mu[i]=rand(-10000..10000)(),i=1..k),seq(alpha[i]=rand(-10000..10000)(),i=1..(k-1))],Jac):
			lower_bound_of_dimension := LinearAlgebra[Rank](%);
		end if:

		if lower_bound_of_dimension <> expected_dimension then
			subs([seq(lambda[i]=rand(1..10000)(),i=1..k),seq(mu[i]=rand(-10000..10000)(),i=1..k),seq(alpha[i]=rand(-10000..10000)(),i=1..(k-1))],Jac):
			lower_bound_of_dimension := LinearAlgebra[Rank](%);
		end if:

		print(lower_bound_of_dimension,expected_dimension);

		if lower_bound_of_dimension <> expected_dimension then
			print("Lower than expected rank")
		end if:

		print(\n);

		break if lower_bound_of_dimension = d;

	end do;
end do;
