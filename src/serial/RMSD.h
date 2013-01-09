#ifndef RMSD_H_
#define RMSD_H_

#include <vector>

class RMSD{

	public:
		RMSD(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		virtual ~RMSD(){};
		virtual void oneVsFollowing(int conformation, double* rmsd) = 0;
		virtual void calculateRMSDCondensedMatrix(std::vector<double>& rmsd) = 0;
		virtual void iterativeSuperposition(double);
	protected:
		RMSD(){}
		int numberOfConformations;
		int atomsPerConformation;
		int coordinatesPerConformation;
		double* allCoordinates;

	private:
		RMSD(const RMSD& r){}
};
#endif
