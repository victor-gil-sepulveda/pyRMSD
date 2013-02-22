# Conditional import
try:
    import prody
    from prody.proteins import parsePDB
    prody_available = True

except ImportError:
    prody_available = False
    
import numpy
import pyRMSD.pdbReader

def flattenCoords(coordsets):
    """
    Flattens a coordinate set in Prody's array format into a 1D array.
    
    @param coordsets: The coodinates in Prody's array format
     
    @author: vgil
    @date: 30/11/2012
    """
    return numpy.reshape(coordsets,coordsets.shape[0]*coordsets.shape[1]*coordsets.shape[2])
    
class Reader(object):
    
    @staticmethod
    def availableReaderTypes(): return ["PRODY_READER","LITE_READER"]
    
    def __init__ (self, readerType):
        """
        Class constructor. Initializes the small DSL props.
        
        @param readerType: One of the types in 'availableReaderTypes'
         
        @author: vgil
        @date: 30/11/2012
        """
        self.readerType = readerType
        self.numberOfFrames = 0
        self.numberOfAtoms = None
        self.onlyCA = False
        self.filesToRead = []
        self.selectionString = ""
        if(not readerType in Reader.availableReaderTypes()):
            print "The reader type ",readerType, " is not an available reader."
            raise KeyError
        if(self.readerType == "PRODY_READER" and not prody_available):
            print "PRODY_READER is not available, changing to  LITE_READER."
            self.readerType = "LITE_READER"
        
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
    
    #azucaaarr
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
    
    def usingThisSelection(self,sel_string):
        """
        Applies a selection to the read. Only used by the "PRODY_READER" type (prody has a great selection engine).
        
        @param sel_string: The path of the file.
         
        @return: itself, to chain calls.
        
        @author: vgil
        @date: 30/11/2012
        """
        if self.readerType == "PRODY_READER":
            self.selectionString = sel_string
        else:
            print "Your current Reader (",self.readerType,") does not allow the use of selections. Please use the prody (http://www.csb.pitt.edu/prody/) reader for this."
    
    def read(self, verbose = False):
        """
        Does the actual reading using the parameters of the DSL.
        
        @param verbose: If switched on some output is generated.
         
        @return: The read coordinates in the way prody generates.
        
        @author: vgil
        @date: 30/11/2012
        """
        coordinates = numpy.array([])
        
        if self.readerType == "PRODY_READER":
            try:
                # Some prody versions doesn't accept this
                prody.setVerbosity('none')
            except:
                pass
            for pdb_path in self.filesToRead:
                if self.onlyCA:
                    pdb = parsePDB(pdb_path, subset='calpha')
                else:
                    pdb = parsePDB(pdb_path)
                
                if verbose: print "PDB parsed (",pdb_path,")"
                
                if self.selectionString != "":
                    selection = pdb.select(self.selectionString)
                    pdb = selection.copy()
                
                coordsets = pdb.getCoordsets()
                self.numberOfFrames += len(coordsets)
                atoms_in_this_trajectory = len(pdb)
                
                if not self.__checkAndUpdateNumberOfAtoms(atoms_in_this_trajectory,pdb_path):
                    return None
                
                if len(coordinates) == 0:
                    coordinates = coordsets
                else:
                    coordinates = numpy.append(coordinates, coordsets, axis=0)
            #coordinates.shape = (self.numberOfFrames,self.numberOfAtoms,3)
            return coordinates
            
        elif self.readerType == "LITE_READER":
            for path in self.filesToRead:
                if self.onlyCA:
                    coordsets, num_frames, atoms_in_this_trajectory = pyRMSD.pdbReader.readPDB(path," CA ")
                else:
                    coordsets, num_frames, atoms_in_this_trajectory = pyRMSD.pdbReader.readPDB(path)
                
                self.numberOfFrames += num_frames
                
                if not self.__checkAndUpdateNumberOfAtoms(atoms_in_this_trajectory,path):
                    return None
                
                if verbose: print "PDB parsed (",path,")"
                
                coordinates = numpy.append(coordinates, coordsets, axis=0)
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
                print "[Error :: Reader.read] The number of atoms of ", path, " is different from the number of atoms of the first trajectory."
                return False
            else:
                return True
        else:
            self.numberOfAtoms = atoms_in_this_trajectory
            return True
    