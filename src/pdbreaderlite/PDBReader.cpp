///////////////////////////////////////////////////////////
/// PDBReader.cpp
///
/// Implementation of the PDBReader class
///
/// \author vgil
/// \author icabeza
/// \author mrivero
/// \date 28/03/2010
///////////////////////////////////////////////////////////

#include "PDBReader.h"

#include <fstream>
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <cmath>
using namespace std;

PDBReader::PDBReader():number_of_models(0),number_of_atoms(0){}

PDBReader::~PDBReader(){}

void PDBReader::check_atoms(int expected, int now){
	if(expected!=now){
		cout<<"Number of atoms in model "<<number_of_models<<" differ from those of model 0 (Expected "<<expected<<", found "<<now<<")"<<endl;
	}
}


///////////////////////////////////////////////////////////////
/// \remarks
/// This function reads a PDB file.
///
/// \param path [In] Path of the PDB file
///
/// \author vgil
/// \date 10/02/2011
///////////////////////////////////////////////////////////////
void PDBReader::read(const char* path, const char* name_filter_c){
	string line;
	ifstream stream;
	number_of_models=  number_of_atoms = 0;
	all_coordinates.clear();

	int atoms_counter = 0;
	this->all_coordinates.clear();

	// Then let's read the file (the file is readable at least, as it used getTotalModels)
	time_t old_time = time(0);

	stream.open(path);
	if (!stream) {
		cout<<"[Error :: read] Impossible to read file: "<<path<<flush<<endl;
	}

	int line_number = 0;
	string name_filter(name_filter_c);

	// Process the file until the end is reached
	while(!stream.eof()){

		// Read the line and put it in 'line'
		getline(stream,line);
		++line_number;

		// Line processing: Record name goes from character 0 to 5 (included)
		string record = line.substr(0,6);

		// Different record types of interest
		if(record == "ATOM  " || record == "HETATM"){
			processATOM(line, name_filter, atoms_counter);
		}
		else if(record == "MODEL "){

			 if(number_of_atoms == 0 && number_of_models ==1){
				 number_of_atoms = atoms_counter;
			 }

			 check_atoms(number_of_atoms, atoms_counter);

			 atoms_counter = 0;

			 processMODEL();
		}
	}

	check_atoms(number_of_atoms, atoms_counter);

	stream.close();
	cout<<line_number<<" lines read in "<<long(time(0) - old_time)<<" seconds."<<endl;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function processes lines tagged with the ENDMDL tag.
///
/// \param line [In] A line tagged with the ENDMDL tag
///
/// \author vgil
/// \date 10/02/2011
///////////////////////////////////////////////////////////////
void PDBReader::processMODEL(){
     number_of_models++;
}

///////////////////////////////////////////////////////////////
/// \remarks
/// This function processes lines tagged with the ATOM tag.
///
/// \param line [In] A line tagged with the ATOM tag
///
/// \author vgil
/// \date 10/02/2011
///////////////////////////////////////////////////////////////
void PDBReader::processATOM(string& line, string& name_filter, int& atoms_counter){
		 string name = line.substr(12,4);
		 if (name_filter == "" || (name_filter!= "" && name == name_filter)){
			 double x = fromString<double>(line.substr(30,8));
			 all_coordinates.push_back(x);
			 double y = fromString<double>(line.substr(38,8));
			 all_coordinates.push_back(y);
			 double z = fromString<double>(line.substr(46,8));
			 all_coordinates.push_back(z);
			 atoms_counter++;
		 }
}
