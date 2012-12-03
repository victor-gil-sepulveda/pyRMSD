#pragma once
#ifndef __L_PDB__
#define __L_PDB__

#include <vector>
#include <iosfwd>
#include <sstream>


class PDBReader{

	public:

				PDBReader  ();
				~PDBReader ();

		void read(const char* path,const char* name_filter_c);

		std::vector<double> all_coordinates;

		int number_of_models;
		int number_of_atoms;
	private:

		PDBReader(const PDBReader&);
		PDBReader& operator=(const PDBReader&);

		void processLine		(std::string& line);
		void processATOM		(std::string& line, std::string& name_filter, int& atoms_counter);
		void processMODEL 	();

		void check_atoms(int expected, int now);

		template<typename T>
		T fromString(const std::string & s){
		   std::istringstream i(s);
		   T x;
		   i >> x;
		   return x;
		}
};
#endif
