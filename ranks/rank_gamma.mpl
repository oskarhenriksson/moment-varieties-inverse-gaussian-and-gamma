# Fix seed for reproducibility
randomize(12345);

random_parameters := rand(1..1000);

for d from 1 to 100 do

	for k from 2 to d do

		printf("d=%d\n",d);
		printf("k=%d\n",k);

		expected_dimension := min([d,2*k + k - 1]);

		parametrization := subs([alpha[k]=1-add(alpha[i],i=1..(k-1))],[
			seq(add(alpha[i]*theta[i]^r*mul((kappa[i]+j),j=1..(r-1)),i=1..k),r=1..d)]):

		Jac := VectorCalculus[Jacobian](parametrization,[seq(alpha[i],i=1..k-1),seq(kappa[i],i=1..k),seq(theta[i],i=1..k)]):

		# Compute the rank at a random integer point (repeat if expected dimension is not attained)
		lower_bound_of_dimension := -1;
		for i from 1 to 10 do
			if lower_bound_of_dimension < expected_dimension then
				theta_star := [seq(alpha[i]=random_parameters(),i=1..k-1),seq(kappa[i]=random_parameters(),i=1..k),seq(theta[i]=random_parameters(),i=1..k)];
				Jac_star := subs(theta_star,Jac):
				lower_bound_of_dimension := LinearAlgebra[Rank](Jac_star);
			end if;
		end do;
		

		if expected_dimension <> lower_bound_of_dimension then
			printf(">> Lower than expected rank! <<");
		else

			printf(convert(theta_star,string));
			printf("\n");
			printf("\n");

		end if;

		break if lower_bound_of_dimension = d;
	end do;


end do;
