import numpy
import sysconfig
import os.path
from build_utils import compile_a_file_collection, Link
import optparse
import collections


##############################################
#####                                  #######
#####     Compiling Options            #######
#####                                  #######
##############################################
CUDA_BASE = "/usr/local/cuda-4.2"           # Base dir for CUDA installation
CUDA_INCLUDE_FOLDER = CUDA_BASE+"/include"         # CUDA headers path
CUDA_LIBRARIES_FOLDER = CUDA_BASE+"/lib64"            # CUDA libs path ( /lib if you're running it in a 32b machine)
CUDA_ARCHITECHTURE = "sm_11"                         # CUDA architecture of your card.
CUDA_LIBRARY = "cudart"
PYTHON_EXTENSION_OPTIONS = "-pthread -fno-strict-aliasing -fmessage-length=0 -O3 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables -fPIC"
PYTHON_INCLUDE_FOLDER = os.path.dirname(sysconfig.get_paths()['include'])
PYTHON_LIBRARY_FOLDER = os.path.dirname(sysconfig.get_paths()['stdlib'])
PYTHON_LIBRARY = "python2.7"
OPENMP_OPTION = "-fopenmp"
###############################################


##############################################
#####                                  #######
#####     Compiling Options            #######
#####           for the Brave          #######
##############################################
NUMPY_INCLUDE =  numpy.get_include()
PYTHON_EXTENSION_LINKING_OPTIONS = "-pthread -shared"
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
    files_to_compile_with_nvcc = {"src/calculators/QCP":["ThRMSDCuda.cu","kernel_functions_cuda.cu"],
                                  "src/calculators/QCP":["ctest.cpp"]}
    
    
    files_to_compile_with_gcc = {"src/calculators/QCP":["ThRMSDSerial.cpp","kernel_functions_serial.cpp","ThRMSDSerialOmp.cpp"],
                                 "src/calculators":["RMSD.cpp"],
                                 "src/calculators/QTRFIT":["RMSDSerial.cpp"],
                                 "src/matrix":["Matrix.cpp","Statistics.cpp"],
                                 "src/python":["pyRMSD.cpp","readerLite.cpp"],
                                 "src/pdbreaderlite":["PDBReader.cpp"]}
    
    files_to_compile_with_gcc_and_openmp = {"src/calculators/QCP":["kernel_functions_omp.cpp"],
                                            "src/calculators/QTRFIT":["RMSDomp.cpp"],
                                            "src/calculators":["RMSDTools.cpp"],
                                            "src/calculators/test":["TestRMSDTools.cpp"]}
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
                    to_produce("condensedMatrix.so")
    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')            
    os.system( linkDSL.getLinkingCommand())
    
    linkDSL = Link().\
                    using("g++").\
                    with_options([PYTHON_EXTENSION_LINKING_OPTIONS]).\
                    using_libs([PYTHON_LIBRARY]).\
                    using_lib_locations([PYTHON_LIBRARY_FOLDER]).\
                    this_object_files([files_to_link["PDBReader"],files_to_link["readerLite"]]).\
                    to_produce("pdbReader.so")
    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')            
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
                        to_produce("calculators.so")
    else:
        linkDSL = Link().\
                        using("g++").\
                        with_options([PYTHON_EXTENSION_LINKING_OPTIONS,OPENMP_OPTION]).\
                        using_libs([PYTHON_LIBRARY]).\
                        using_lib_locations([PYTHON_LIBRARY_FOLDER]).\
                        this_object_files([files_to_link["ThRMSDSerial"],files_to_link["ThRMSDSerialOmp"],files_to_link["kernel_functions_serial"],\
                                           files_to_link["kernel_functions_omp"],files_to_link["pyRMSD"],\
                                           files_to_link["RMSDomp"],files_to_link["RMSD"],files_to_link["RMSDSerial"],files_to_link["RMSDTools"],]).\
                        to_produce("calculators.so")
    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
    os.system(linkDSL.getLinkingCommand())
    
#    if options.use_cuda:
#        linkDSL = Link().\
#                        using("g++").\
#                        with_options([]).\
#                        using_libs([CUDA_LIBRARY]).\
#                        using_lib_locations([CUDA_LIBRARIES_FOLDER]).\
#                        this_object_files([files_to_link["ctest"],files_to_link["ThRMSDSerial"],files_to_link["kernel_functions_serial"],files_to_link["ThRMSDCuda"],\
#                                           files_to_link["kernel_functions_cuda"],files_to_link["RMSD"]]).\
#                        to_produce("test_main")
#                        
#        os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
#        os.system(linkDSL.getLinkingCommand())

    linkDSL = Link().\
                    using("g++").\
                    with_options([OPENMP_OPTION]).\
                    using_libs([]).\
                    using_lib_locations([]).\
                    this_object_files([files_to_link["TestRMSDTools"],files_to_link["RMSDTools"],files_to_link["RMSDomp"],
                                       files_to_link["RMSD"]]).\
                    to_produce("test_rmsdtools_main")
                    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
    os.system(linkDSL.getLinkingCommand())
    
    os.system("mv calculators.so pyRMSD/")
    os.system("mv condensedMatrix.so pyRMSD/")
    os.system("mv pdbReader.so pyRMSD/")
    os.system("mv test_rmsdtools_main src/calculators/test")
    
#    if options.use_cuda:
#        os.system("mv test_main src/theobald/test/")
    
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
    os.system('echo "\033[33mWriting available calculators...\033[0m"')
    open("pyRMSD/availableCalculators.py","w").write(calcs_str)
    
    os.system('echo "\033[32mCleaning...\033[0m"')
    for produced_file in  files_to_link:
        if files_to_link[produced_file] != "":
            os.system("rm "+files_to_link[produced_file])

