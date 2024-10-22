# Fix seed for reproducibility
randomize(12345);

random_parameters := rand(1..1000);

for d from 1 to 10 do
	for k from 2 to d do

		printf("d=%d\n",d);
		printf("k=%d\n",k);

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
		# Make 10 attempts to obtain the expected dimension

		lower_bound_of_dimension := -1;

		for i from 1 to 10 do			
			if lower_bound_of_dimension <> expected_dimension then
				theta_star := [seq(alpha[i]=random_parameters(),i=1..(k-1)),seq(lambda[i]=random_parameters(),i=1..k),seq(mu[i]=random_parameters(),i=1..k)]:
				Jac_star := subs(theta_star,Jac):
				lower_bound_of_dimension := LinearAlgebra[Rank](Jac_star);
			end if:
		end do;


		#print(lower_bound_of_dimension,expected_dimension);

		if lower_bound_of_dimension <> expected_dimension then
			printf(">> Lower than expected rank! <<")
		else
			printf(convert(theta_star,string));
			printf("\n");
			printf("\n");
		end if;

		break if lower_bound_of_dimension = d;

	end do;
end do;
