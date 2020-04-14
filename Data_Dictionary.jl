using DelimitedFiles

function DataDictionary(time_start,time_stop,time_step)

	# Load the stoichiometric network from disk -
	stoichiometric_matrix = readdlm("C:\\Users\\user\\Desktop\\Cornell Academics\\Spring 2020\\CHEME 7770\\PS3\\S_Matrix.dat");

	#calculate vmax values
	kcat_v1 = 203 * 3600;		#converting from s^-1 to hr^-1
	kcat_v2 = 34.5 * 3600;
	kcat_v3 = 249 * 3600;
	kcat_v4 = 88.1 * 3600;
	kcat_v5 = 13.7 * 3600;

	E = 0.01*10^-3;		#converting from umol/gdw to mmol/gdw

	vmax_1 = kcat_v1 * E;
	vmax_2 = kcat_v2 * E;
	vmax_3 = kcat_v3 * E;
	vmax_4 = kcat_v4 * E;
	vmax_5 = kcat_v5 * E;

	# Saturation function: [S]/(Km+[S]) when available
	sat_1 = (4.67E-3/(4.67E-3 + 3.92E-4))*(1.49E-2/(1.49E-2 + 1.54E-4));
	sat_2 = 1;
	sat_3 = 2.55E-4/(1.55E-3 + 2.55E-4)	;
	sat_4 = 4.49E-3/(1.6E-3 + 4.49E-3);
	sat_5 = 2.55E-4/(4.4E-03 + 2.55E-4);

	# Setup default flux bounds array -
	default_bounds_array = [

			0	vmax_1*sat_1;		# v1	ATP + L-citrulline + L-aspartate --> AMP + diphosphate + 2-(Nomega-L-arginino)succinate
			0	vmax_2*sat_2;		# v2 	2-(Nomega-L-arginino)succinate --> fumarate + L-arginine
			0	vmax_3*sat_3; 		# v3 	L-arginine + H2O --> L-ornithine + urea
			0	vmax_4*sat_4;		# v4 	carbamoyl phosphate + L-ornithine --> phosphate + L-citrulline
 -vmax_5*sat_5	vmax_5*sat_5;		# v5 	2 L-arginine + 3 NADPH + 3 H+ + 4 O2 <--> 2 L-citrulline + 2 nitric oxide + 3 NADP+ + 4 H2O
			0	10;		# b1 	[] --> carbomyl phosphate
			0	10;		# b2 	[] --> aspartate
			0	10;		# b3  	fumarate --> []
			0	10;		# b4 	urea --> []
			0	10;		# b5 	[] --> ATP
			0	10;		# b6 	AMP --> []
			0	10;		# b7 	diphosphate --> []
		    0	10;		# b8 	[] --> H2O
			0	10;		# b9 	phosphate --> []
		  -10	10;		# b10 	[] --> NADPH --> []
		  -10	10;		# b11 	[] --> H+ --> []
		  -10	10;		# b12 	[] --> O2 --> []
		  -10	10;		# b13 	[] --> NO --> []
		  -10	10;		# b14 	[] --> NADP+ --> []
			0	10;		# b15 	H2O --> []
		];

	# Setup default species bounds array -
	species_bounds_array = [

		0.0	0.0	;	# 1 aspartate
		0.0	0.0	;	# 2 ATP
		0.0	0.0	;	# 3 citrulline
		0.0	0.0	;	# 4 AMP
		0.0	0.0	;	# 5 diphosphate
		0.0	0.0	;	# 6 argininosuccinate
		0.0	0.0	;	# 7 fumarate
		0.0	0.0	;	# 8 arginine
		0.0	0.0	;	# 9 H2O
		0.0	0.0	;	# 10 ornithine
		0.0	0.0	;	# 11 urea
		0.0	0.0	;	# 12 carbomyl phosphate
		0.0	0.0	;	# 13 phosphate
		0.0	0.0	;	# 14 NADPH
		0.0	0.0	;	# 15 H+
		0.0	0.0	;	# 16 O2
		0.0	0.0	;	# 17 NO
		0.0	0.0	;	# 18 NADP+

	];

	# Setup the objective coefficient array -
	objective_coefficient_array = [

		0	; # v1
		0	; # v2
		0	; # v3
		0	; # v4
		0	; # v5
		0	; # b1
		0	; # b2
		0	; # b3
		1	; # b4
		0	; # b5
		0	; # b6
		0	; # b7
		0	; # b8
		0	; # b9
		0	; # b10
		0	; # b11
		0	; # b12
		0	; # b13
		0	; # b14
		0	; # b15
	];

	# Min/Max flag - default is minimum -
	min_flag = false

	# List of reation strings - used to write flux report
	list_of_reaction_strings = [

	"v1: ATP + L-citrulline + L-aspartate --> AMP + diphosphate + 2-(Nomega-L-arginino)succinate"
	"v2: 2-(Nomega-L-arginino)succinate --> fumarate + L-arginine"
	"v3: L-arginine + H2O --> L-ornithine + urea"
	"v4: carbamoyl phosphate + L-ornithine --> phosphate + L-citrulline"
	"v5: 2 L-arginine + 3 NADPH + 3 H+ + 4 O2 --> 2 L-citrulline + 2 nitric oxide + 3 NADP+ + 4 H2O"
	"b1: [] --> carbomyl phosphate"
	"b2: [] --> aspartate"
	"b3: fumarate --> []"
	"b4: urea --> []"
	"b5: [] --> ATP"
	"b6: AMP --> []"
	"b7: diphosphate --> []"
	"b8: [] --> H2O"
	"b9: phosphate --> []"
	"b10: [] --> NADPH --> []"
	"b11: [] --> H+ --> []"
	"b12: [] --> O2 --> []"
	"b13: [] --> NO --> []"
	"b14: [] --> NADP+ --> []"
	"b15: H2O --> []"

	];

	# List of metabolite strings - used to write flux report
	list_of_metabolite_symbols = [

		"aspartate"
		"ATP"
		"citrulline"
		"AMP"
		"diphosphate"
		"argininosuccinate"
		"fumarate"
		"arginine"
		"H2O"
		"ornithine"
		"urea"
		"carbamoyl phosphate"
		"phosphate"
		"NADPH"
		"H+"
		"O2"
		"NO"
		"NADP+"

	];

	# =============================== DO NOT EDIT BELOW THIS LINE ============================== #
	data_dictionary = Dict{AbstractString,Any}()
	data_dictionary["stoichiometric_matrix"] = stoichiometric_matrix
	data_dictionary["objective_coefficient_array"] = objective_coefficient_array
	data_dictionary["default_flux_bounds_array"] = default_bounds_array;
	data_dictionary["species_bounds_array"] = species_bounds_array
	data_dictionary["list_of_reaction_strings"] = list_of_reaction_strings
	data_dictionary["list_of_metabolite_symbols"] = list_of_metabolite_symbols
	data_dictionary["min_flag"] = min_flag
	# =============================== DO NOT EDIT ABOVE THIS LINE ============================== #
	return data_dictionary
end
