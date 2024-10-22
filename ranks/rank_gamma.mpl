# Fix seed for reproducibility
randomize(12345);

for d from 1 to 100 do
	for k from 1 to d do

		print((d,k));

		expected_dimension := min([d,2*k + k - 1]);

		parametrization := subs([lambda[k]=1-add(lambda[i],i=1..(k-1))],[
			seq(add(lambda[i]*beta[i]^r*mul((alpha[i]+j),j=1..(r-1)),i=1..k),r=1..d)]):

		Jac := VectorCalculus[Jacobian](parametrization,[seq(lambda[i],i=1..k-1),seq(alpha[i],i=1..k),seq(beta[i],i=1..k)]):

		# Compute the rank at a random integer point
		subs([seq(lambda[i]=rand(-1000..1000)(),i=1..k-1),seq(alpha[i]=rand(-1000..1000)(),i=1..k),seq(beta[i]=rand(-1000..1000)(),i=1..k)],Jac):
		lower_bound_of_dimension := LinearAlgebra[Rank](%);

		# Pick another point if the rank is lower than the expected one
		if expected_dimension <> lower_bound_of_dimension then
			subs([seq(lambda[i]=rand(-10000..10000)(),i=1..k-1),seq(alpha[i]=rand(-10000..10000)(),i=1..k),seq(beta[i]=rand(-10000..10000)(),i=1..k)],Jac):
			lower_bound_of_dimension := LinearAlgebra[Rank](%);
		end if;

		if expected_dimension <> lower_bound_of_dimension then
			subs([seq(lambda[i]=rand(-10000..10000)(),i=1..k-1),seq(alpha[i]=rand(-10000..10000)(),i=1..k),seq(beta[i]=rand(-10000..10000)(),i=1..k)],Jac):
			lower_bound_of_dimension := LinearAlgebra[Rank](%);
		end if;

		print((lower_bound_of_dimension,expected_dimension));

		if expected_dimension <> lower_bound_of_dimension then
			print("Lower than expected rank");
		end if;

		print(\n);

		break if lower_bound_of_dimension = d;
	end do;
end do;
