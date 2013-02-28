#ifndef RMSD_H_
#define RMSD_H_

#include <vector>

class RMSD{

	public:
		RMSD(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		void setCalculationCoordinates(int atomsPerRMSDConformation, double* const allRMSDCoordinates);
		virtual ~RMSD(){};
		virtual void oneVsFollowing(int conformation, double* rmsd) = 0;
		virtual void calculateRMSDCondensedMatrix(std::vector<double>& rmsd) = 0;
		virtual void iterativeSuperposition(double);

	protected:
		int numberOfConformations;

		int atomsPerConformation;
		int coordinatesPerConformation;
		double* allCoordinates; // Coordinates for fitting and RMSD (if allRMSDCoordinates == NULL)

		int atomsPerRMSDConformation;
		int coordinatesPerRMSDConformation;
		double* allRMSDCoordinates; // If is different from NULL, then this are the coordinates
									 // to calculate the RMSD.

		bool modifyFittingCoordinates;

	private:
		RMSD(){}
		RMSD(const RMSD& r){}
};
#endif
