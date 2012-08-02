#ifndef RMSD_H_
#define RMSD_H_

#include <vector>

class RMSD{

	public:
		RMSD(int numberOfConformations, int atomsPerConformation, double* allCoordinates);
		virtual ~RMSD(){};
		virtual void oneVsTheOthers(int conformation, double* rmsd) = 0;
		virtual void calculateRMSDCondensedMatrix(std::vector<double>& rmsd) = 0;

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
