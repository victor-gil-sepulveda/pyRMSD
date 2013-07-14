#include "tests.h"
#include "../factory/RMSDCalculatorTypes.h"
#include <vector>
#include <cstdlib>
using std::vector;


void set_precision_if(
		RMSDCalculatorType calctype,
		double alternative_precision,
		double else_precision,
		double& precision){
	precision = else_precision;
	#ifdef USE_CUDA
		#ifdef CUDA_PRECISION_SINGLE
		if (calctype == QCP_CUDA_CALCULATOR)
			precision = alternative_precision;
		#endif
	#endif
}

int main(int argc, char **argv){

	RMSDCalculatorType available_calculators_d [] =  {
			KABSCH_SERIAL_CALCULATOR,
//			KABSCH_OMP_CALCULATOR,
//			QTRFIT_SERIAL_CALCULATOR,
//			QTRFIT_OMP_CALCULATOR,
//			QCP_SERIAL_CALCULATOR,
//			QCP_OMP_CALCULATOR,
#ifdef USE_CUDA
			QCP_CUDA_CALCULATOR
#endif
	};

	vector<RMSDCalculatorType> available_calculators( available_calculators_d,
			available_calculators_d + sizeof(available_calculators_d)/sizeof(RMSDCalculatorType));

	test_initialize();

	test_copy_array();

	test_mean_coordinates();

	test_translations();

	test_center_coordinates();

	test_QCP_Kernel();

	test_KABSCH_Kernel();

	for(unsigned int i = 0; i < available_calculators.size();++i){

		test_superposition_with_fit(available_calculators[i],
				"test_data/Models/prot_plus_ligand_similar/prot_plus_ligand_offset.CA.coords",
				"test_data/Superpose_Fit_CA_similar/prot_plus_ligand_similar.aligned_CA.coords",
				"test_data/Superpose_Fit_CA_similar/prot_plus_ligand_similar.aligned_CA.rmsd",
				1e-12);

		test_superposition_with_fit(available_calculators[i],
				"test_data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords",
				"test_data/Superpose_Fit_CA_very_diff/prot_plus_ligand.aligned_CA.coords",
				"test_data/Superpose_Fit_CA_very_diff/prot_plus_ligand.aligned_CA.rmsd",
				1e-12);

		test_superposition_with_fit_and_calc(available_calculators[i],
				"test_data/Models/prot_plus_ligand_similar/prot_plus_ligand_offset.CA.coords",
				"test_data/Superpose_Fit_CA_Calc_BEN_similar/prot_plus_ligand_similar.aligned_CA.coords",
				"test_data/Models/prot_plus_ligand_similar/prot_plus_ligand_offset.ligand.coords",
				"test_data/Superpose_Fit_CA_Calc_BEN_similar/prot_plus_ligand_similar.aligned_BEN.coords",
				"test_data/Superpose_Fit_CA_Calc_BEN_similar/prot_plus_ligand_similar.aligned_BEN.rmsd",
				1e-12);

		test_superposition_with_fit_and_calc(available_calculators[i],
				"test_data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.CA.coords",
				"test_data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_CA.coords",
				"test_data/Models/prot_plus_ligand_very_different/not_aligned_offset_prot_plus_ligand.ligand.coords",
				"test_data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.coords",
				"test_data/Superpose_Fit_CA_Calc_BEN_very_different/prot_plus_ligand_similar.aligned_BEN.rmsd",
				1e-12);

		test_step_by_step_iterative_superposition_with_fit(available_calculators[i],
						"test_data/Iterpose_Fit_CA/steps",
						"test_data/Iterpose_Fit_CA/mean",
						"test_data/Iterpose_Fit_CA/stretching_trajectory_offset_ligand.initial_CA.coords",
						1e-12,
						9);

		test_iterative_superposition_with_fit(available_calculators[i],
						"test_data/Iterpose_Fit_CA/stretching_trajectory_offset_ligand.initial_CA.coords",
						"test_data/Iterpose_Fit_CA/stretching_trajectory_offset_ligand.iterposed_CA.coords",//steps/iter_step_8.coords",
						"test_data/Iterpose_Fit_CA/step_rmsd_diff.rmsd",
						1e-12,
						9);

		test_iterative_superposition_with_fit_and_calc_rotation(available_calculators[i],
						"test_data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.initial_all.coords",
						"test_data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.iterposed_all.coords",
						"test_data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.initial_BEN.coords",
						"test_data/Iterpose_Fit_CA_Rot_BEN/stretching_trajectory_offset_ligand.iterposed_BEN.coords",
						"test_data/Iterpose_Fit_CA_Rot_BEN/step_rmsd_diff.rmsd",
						1e-12,
						10);
	}

	return 0;
}

