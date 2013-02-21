import numpy
import sysconfig
import os.path
from build_utils import compile_a_file_collection, Link
import optparse
import collections
import shutil

##############################################
#####                                  #######
#####     Compiling Options            #######
#####                                  #######
##############################################
CUDA_BASE = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.0"#"/usr/local/cuda-4.2"                   # Base dir for CUDA installation
CUDA_INCLUDE_FOLDER = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.0\include" #CUDA_BASE+"/include"          # CUDA headers path
CUDA_LIBRARIES_FOLDER = "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v5.0\lib\x64" #CUDA_BASE+"/lib64"          # CUDA libs path ( /lib if you're running it in a 32b machine)
CUDA_ARCHITECHTURE = "sm_11"                        # CUDA architecture of your card.
CUDA_LIBRARY = "cudart"
PYTHON_EXTENSION_OPTIONS = "-mdll -O -Wall -fno-strict-aliasing -fmessage-length=0 -O3 -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables"
PYTHON_INCLUDE_FOLDER = "C:\Python27\include" #os.path.dirname(sysconfig.get_paths()['include'])
PYTHON_LIBRARY_FOLDER = "C:\Python27\libs" #os.path.dirname(sysconfig.get_paths()['stdlib'])
PYTHON_LIBRARY = "python27"
OPENMP_OPTION = "-fopenmp"
###############################################


##############################################
#####                                  #######
#####     Compiling Options            #######
#####           for the Brave          #######
##############################################
NUMPY_INCLUDE =  numpy.get_include()
PYTHON_EXTENSION_LINKING_OPTIONS = "-shared -s "
CUDA_OPTIONS = "-O3 -use_fast_math --gpu-architecture %s --compiler-options '-fPIC'"%(CUDA_ARCHITECHTURE)
DEFINE_USE_CUDA = "-DUSE_CUDA"
BASE_DIR = os.getcwd()
###############################################

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog [--cuda]', version='1.0')
    parser.add_option('--cuda',  action="store_true",dest = "use_cuda",  help="Use this flag if you want to compile the CUDA calculator.")
    options, args = parser.parse_args()
    
    ######################################
    ### Files we are going to use
    ######################################
    files_to_compile_with_nvcc = {"src\\theobald":["ThRMSDCuda.cu","kernel_functions_cuda.cu"],
                                  "src\\theobald\\test":["ctest.cpp"]}
    
    
    files_to_compile_with_gcc = {"src\\theobald":["ThRMSDSerial.cpp","kernel_functions_serial.cpp","ThRMSDSerialOmp.cpp"],
                                 "src\\serial":["RMSDSerial.cpp","RMSD.cpp"],
                                 "src\\matrix":["Matrix.cpp","Statistics.cpp"],
                                 "src\\python":["pyRMSD.cpp","readerLite.cpp"],
                                 "src\\pdbreaderlite":["PDBReader.cpp"]}
    
    files_to_compile_with_gcc_and_openmp = {"src\\theobald":["kernel_functions_omp.cpp"],
                                            "src\\serial":["RMSDomp.cpp","RMSDTools.cpp"],
                                            "src\\serial\\test":["TestRMSDTools.cpp"]}
    #########################################
    
    files_to_link = collections.defaultdict(str)
    
    if options.use_cuda:
        compile_a_file_collection(BASE_DIR, files_to_compile_with_gcc, "gcc", PYTHON_EXTENSION_OPTIONS+" "+DEFINE_USE_CUDA , [PYTHON_INCLUDE_FOLDER, NUMPY_INCLUDE], ".o",files_to_link)
    else:
        compile_a_file_collection(BASE_DIR, files_to_compile_with_gcc, "gcc", PYTHON_EXTENSION_OPTIONS, [PYTHON_INCLUDE_FOLDER, NUMPY_INCLUDE], ".o",files_to_link)
        
    compile_a_file_collection(BASE_DIR, files_to_compile_with_gcc_and_openmp, "gcc", PYTHON_EXTENSION_OPTIONS+" "+OPENMP_OPTION, [PYTHON_INCLUDE_FOLDER, NUMPY_INCLUDE], ".o",files_to_link)
    
    if options.use_cuda:
        compile_a_file_collection(BASE_DIR, files_to_compile_with_nvcc, "nvcc", CUDA_OPTIONS, [CUDA_INCLUDE_FOLDER], ".o",files_to_link)
    
    linkDSL = Link().\
                    using("g++").\
                    with_options([PYTHON_EXTENSION_LINKING_OPTIONS,OPENMP_OPTION]).\
                    using_libs([PYTHON_LIBRARY]).\
                    using_lib_locations([PYTHON_LIBRARY_FOLDER]).\
                    this_object_files([files_to_link["Matrix"],files_to_link["Statistics"]]).\
                    to_produce("condensedMatrix.pyd")
    
    os.system('echo '+ linkDSL.getLinkingCommand()+'')            
    os.system( linkDSL.getLinkingCommand())
    
    linkDSL = Link().\
                    using("g++").\
                    with_options([PYTHON_EXTENSION_LINKING_OPTIONS]).\
                    using_libs([PYTHON_LIBRARY]).\
                    using_lib_locations([PYTHON_LIBRARY_FOLDER]).\
                    this_object_files([files_to_link["PDBReader"],files_to_link["readerLite"]]).\
                    to_produce("pdbReader.pyd")
    
    os.system('echo '+ linkDSL.getLinkingCommand())            
    os.system( linkDSL.getLinkingCommand())
    
    if options.use_cuda:
        linkDSL = Link().\
                        using("g++").\
                        with_options([PYTHON_EXTENSION_LINKING_OPTIONS,OPENMP_OPTION]).\
                        using_libs([PYTHON_LIBRARY,CUDA_LIBRARY]).\
                        using_lib_locations([CUDA_LIBRARIES_FOLDER,PYTHON_LIBRARY_FOLDER]).\
                        this_object_files([files_to_link["ThRMSDSerial"],files_to_link["ThRMSDSerialOmp"],files_to_link["ThRMSDCuda"],files_to_link["kernel_functions_serial"],\
                                           files_to_link["kernel_functions_omp"],files_to_link["kernel_functions_cuda"],files_to_link["pyRMSD"],\
                                           files_to_link["RMSDomp"],files_to_link["RMSD"],files_to_link["RMSDSerial"],files_to_link["RMSDTools"],]).\
                        to_produce("calculators..pyd")
    else:
        linkDSL = Link().\
                        using("g++").\
                        with_options([PYTHON_EXTENSION_LINKING_OPTIONS,OPENMP_OPTION]).\
                        using_libs([PYTHON_LIBRARY]).\
                        using_lib_locations([PYTHON_LIBRARY_FOLDER]).\
                        this_object_files([files_to_link["ThRMSDSerial"],files_to_link["ThRMSDSerialOmp"],files_to_link["kernel_functions_serial"],\
                                           files_to_link["kernel_functions_omp"],files_to_link["pyRMSD"],\
                                           files_to_link["RMSDomp"],files_to_link["RMSD"],files_to_link["RMSDSerial"],files_to_link["RMSDTools"],]).\
                        to_produce("calculators.pyd")
    
    os.system('echo '+ linkDSL.getLinkingCommand())
    os.system(linkDSL.getLinkingCommand())
    
    linkDSL = Link().\
                    using("g++").\
                    with_options([OPENMP_OPTION]).\
                    using_libs([]).\
                    using_lib_locations([]).\
                    this_object_files([files_to_link["TestRMSDTools"],files_to_link["RMSDTools"],files_to_link["RMSDomp"],
                                       files_to_link["RMSD"]]).\
                    to_produce("test_rmsdtools_main.exe")
                    
    os.system('echo '+ linkDSL.getLinkingCommand())
    os.system(linkDSL.getLinkingCommand())
    
    for product in ["pyRMSD\\calculators.pyd","pyRMSD\\condensedMatrix.pyd","pyRMSD\\pdbReader.pyd","src\\serial\\test\\test_rmsdtools_main.exe"]:
        os.system("del "+product)

    shutil.move("calculators.pyd", "pyRMSD")
    shutil.move("condensedMatrix.pyd", "pyRMSD") 
    shutil.move("pdbReader.pyd", "pyRMSD")  
    #shutil.move("test_rmsdtools_main", "src\serial\test")
    
    os.system('echo Cleaning...')
    for produced_file in  files_to_link:
        if files_to_link[produced_file] != "":
            os.remove(files_to_link[produced_file])
            
    ##Calculators
    if options.use_cuda:
        calcs_str = """
def availableCalculators():
    return {"KABSCH_PYTHON_CALCULATOR":-1,"QTRFIT_SERIAL_CALCULATOR":0,"QTRFIT_OMP_CALCULATOR":1,"QCP_CUDA_CALCULATOR":2,"QCP_SERIAL_CALCULATOR":3,"QCP_OMP_CALCULATOR":4}
"""
    else:
        calcs_str = """
def availableCalculators():
    return {"KABSCH_PYTHON_CALCULATOR":-1,"QTRFIT_SERIAL_CALCULATOR":0,"QTRFIT_OMP_CALCULATOR":1,"QCP_SERIAL_CALCULATOR":3,"QCP_OMP_CALCULATOR":4}
"""
    
    os.system('echo Writing available calculators...')
    open("pyRMSD/availableCalculators.py","w").write(calcs_str)
    
