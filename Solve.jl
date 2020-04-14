include("Flux.jl")
include("Data_Dictionary.jl")

data_dictionary = DataDictionary(0,0,0)
#calculate_optimal_flux_distribution(data_dictionary["stoichiometric_matrix"],data_dictionary["default_flux_bounds_array"],data_dictionary["species_bounds_array"],data_dictionary["objective_coefficient_array"], data_dictionary["min_flag"])

(objective_value, calculated_flux_array, dual_value_array, uptake_array, exit_flag, status_flag) = calculate_optimal_flux_distribution(data_dictionary["stoichiometric_matrix"],data_dictionary["default_flux_bounds_array"],data_dictionary["species_bounds_array"],data_dictionary["objective_coefficient_array"], data_dictionary["min_flag"])

print("Calculated flux area : ")
display(calculated_flux_array)
print("Maximum rate of urea production : ",calculated_flux_array[9])
