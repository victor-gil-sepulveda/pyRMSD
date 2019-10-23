from __future__ import print_function
import numpy
import pyRMSD.pdbReader

class Reader(object):

    def __init__ (self):
        """
        Class constructor. Initializes the small DSL props.
        
        @author: vgil
        @date: 30/11/2012
        """
        self.numberOfFrames = 0
        self.numberOfAtoms = None
        self.onlyCA = False
        self.filesToRead = []
        
    def readThisFile(self, file_path):
        """
        Specifies the file to be read.
        
        @param file_path: The path of the file.
         
        @return: itself, to chain calls.
        
        @author: vgil
        @date: 30/11/2012
        """
        self.filesToRead.append(file_path)
        return self
    
    # syntactical sugar :)
    def andThisOtherFile(self, file_path):
        """
        Syntactic sugar.
        @see: readThisFile
        
        @param file_path: The path of the file.
         
        @return: itself, to chain calls.
        
        @author: vgil
        @date: 30/11/2012
        """
        return self.readThisFile(file_path)
    
    def gettingOnlyCAs(self):
        """
        If used, the reader will only read CAs (speeding up the reading process a lot).
        
        @return: itself, to chain calls.
        
        @author: vgil
        @date: 30/11/2012
        """
        self.onlyCA = True
        return self
    
    def read(self, verbose = False):
        """
        Does the actual reading using the parameters of the DSL.
        
        @param verbose: If switched on some output is generated.
         
        @return: The read coordinates in the way prody generates.
        
        @author: vgil
        @date: 30/11/2012
        """
        coordinates = numpy.array([])
        for path in self.filesToRead:
            reader = pyRMSD.pdbReader.PDBReader()
            if self.onlyCA:
                coordsets, num_frames, atoms_in_this_trajectory = reader.read(path," CA ")
            else:
                coordsets, num_frames, atoms_in_this_trajectory = reader.read(path)
            
            self.numberOfFrames += num_frames
            
            if not self.__checkAndUpdateNumberOfAtoms(atoms_in_this_trajectory,path):
                return None
            
            if verbose: print("PDB parsed (",path,")")
            
            coordinates = numpy.append(coordinates, coordsets, axis=0)
            del coordsets
            del reader
        print(coordinates.shape, (self.numberOfFrames, self.numberOfAtoms, 3))
        coordinates.shape = (self.numberOfFrames, self.numberOfAtoms, 3)
        return coordinates
        
    def __checkAndUpdateNumberOfAtoms(self, atoms_in_this_trajectory, path):
        """
        Checks if the new read trajectory has models with the same number of atoms.
        
        @param atoms_in_this_trajectory: The number of atoms of the newly read trajectory.
        
        @return: True or false depending if coherence is maintained or not.
        
        @author: vgil
        @date: 30/11/2012
        """
        if( self.numberOfAtoms != None):
            if self.numberOfAtoms != atoms_in_this_trajectory:
                print("[Error :: Reader.read] The number of atoms of ", path, " is different from the number of atoms of the first trajectory.")
                return False
            else:
                return True
        else:
            self.numberOfAtoms = atoms_in_this_trajectory
            return True
    
