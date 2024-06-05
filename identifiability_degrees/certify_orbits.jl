using DelimitedFiles
using Dates

function certify_orbits(F::System, R::MonodromyResult, action; save_results=false)
    
    representatives = solution.(R.results)
    parameter_values = R.parameters
    
    if save_results
        matrix_of_representatives = hcat(representatives...)
        current_datetime = now()
        iso_datetime_string = Dates.format(current_datetime, "yyyy-mm-ddTHH:MM:SS")
        writedlm("matrix_of_representatives_"*iso_datetime_string*".txt", matrix_of_representatives, ' ')
        writedlm("choice_of_parameters_"*iso_datetime_string*".txt", R.parameters, ' ')
    end
    
    # Certify the representatives
    C = certify(F, representatives, target_parameters=parameter_values)
    certified_representatives = [solution_approximation(c) for c in certificates(C) if is_certified(c)]
    number_of_certified_representatives = length(certified_representatives)
    
    if save_results
        save("representative_certificates_"*iso_datetime_string*".txt",C)
    end
    
    # For each certifiable representative, check that the whole orbit is certifiable
    representatives_of_certifiable_orbits = Vector{ComplexF64}[]
    union_of_certifiable_orbits = Vector{ComplexF64}[]
    for (i,c) in enumerate( certified_representatives )
        if i%100==0
            @info "Currently certifying: orbit $(i) out of $(number_of_certified_representatives)"
        end
        full_orbit = [[v...] for v in action(c)]
        C_c = certify(F, full_orbit, target_parameters=parameter_values)
        if ndistinct_certified(C_c)==length(full_orbit)
            push!(representatives_of_certifiable_orbits,c)
            append!(union_of_certifiable_orbits,full_orbit)
        end
    end

    number_of_certified_orbits = length(representatives_of_certifiable_orbits)

    # Check that all of the solutions are certifiably distinct (so that the oribts don't overlap)
    C_union = certify(F, union_of_certifiable_orbits, target_parameters=parameter_values)
    
    if save_results
        save("union_certificates"*iso_datetime_string,C_union)
    end
    
    if ndistinct_certified(C_union) == length(union_of_certifiable_orbits)
        println("Found $(number_of_certified_orbits) certifiable and nonoverlapping orbits.")
        println("$(number_of_certified_representatives-number_of_certified_orbits) original orbits were discarded")
    else
        error("Overlapping orbits.")
    end

end
