import numpy
import sysconfig
import os
from install_utils import compile_a_file_collection, Link
import optparse
import collections


##############################################
#####                                  #######
#####     Compiling Options            #######
#####                                  #######
##############################################


CUDA_BASE = "/usr/local/cuda-4.2"           # Base dir for CUDA installation
CUDA_INCLUDE = CUDA_BASE+"/include"         # CUDA headers path
CUDA_LIB_64 = CUDA_BASE+"/lib64"            # CUDA libs path ( /lib if you're running it in a 32b machine)
CUDA_ARCH = "sm_11"                         # CUDA architecture of your card.

PYTHON_EXTENSION_OPTIONS = "-pthread -fno-strict-aliasing -fmessage-length=0 -O2 -Wall -D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables -fPIC"
OPENMP_OPTION = "-fopenmp"
PYTHON_LIBRARY = "python2.7"
###############################################


##############################################
#####                                  #######
#####     Compiling Options            #######
#####           for the Brave          #######
##############################################
PYTHON_EXTENSION_LINKING_OPTIONS = "-pthread -shared"
PYTHON_INCLUDE = sysconfig.get_paths()['platinclude']
NUMPY_INCLUDE =  numpy.get_include()
DEFINE_USE_CUDA = "-DUSE_CUDA"
CUDA_OPTIONS = "-O3 -use_fast_math --gpu-architecture %s --compiler-options '-fPIC'"%(CUDA_ARCH)
BASE_DIR = os.getcwd()
###############################################


if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog -m <arg> -c <arglist> [-o <arg>]', version='1.0')
    parser.add_option('--cuda',  action="store_true",dest = "use_cuda",  help="Use this flag if you want to compile the CUDA calculator.")
    options, args = parser.parse_args()
    
    ######################################
    ### Files we are going to use
    ######################################
    files_to_compile_with_nvcc = {"src/theobald":["ThRMSDCuda.cu","kernel_functions_cuda.cu"],
                                  "src/theobald/test":["ctest.cpp"]}
    
    
    files_to_compile_with_gcc = {"src/theobald":["ThRMSDSerial.cpp","kernel_functions_serial.cpp","ThRMSDSerialOmp.cpp"],
                                 "src/serial":["RMSDSerial.cpp","RMSD.cpp","RMSDTools.cpp"],
                                 "src/matrix":["matrix.cpp"],
                                 "src/python":["pyRMSD.cpp","NumpyHelperFuncs.cpp"]}
    
    files_to_compile_with_gcc_and_openmp = {"src/theobald":["kernel_functions_omp.cpp"],
                                            "src/serial":["RMSDomp.cpp"]}
    #########################################
    
    files_to_link = collections.defaultdict(str)
    
    if options.use_cuda:
        compile_a_file_collection(BASE_DIR, files_to_compile_with_gcc, "gcc", PYTHON_EXTENSION_OPTIONS+" "+DEFINE_USE_CUDA , [PYTHON_INCLUDE, NUMPY_INCLUDE], ".o",files_to_link)
    else:
        compile_a_file_collection(BASE_DIR, files_to_compile_with_gcc, "gcc", PYTHON_EXTENSION_OPTIONS, [PYTHON_INCLUDE, NUMPY_INCLUDE], ".o",files_to_link)
        
    compile_a_file_collection(BASE_DIR, files_to_compile_with_gcc_and_openmp, "gcc", PYTHON_EXTENSION_OPTIONS+" "+OPENMP_OPTION, [PYTHON_INCLUDE, NUMPY_INCLUDE], ".o",files_to_link)
    
    if options.use_cuda:
        compile_a_file_collection(BASE_DIR, files_to_compile_with_nvcc, "nvcc", CUDA_OPTIONS, [CUDA_INCLUDE], ".o",files_to_link)
    
    linkDSL = Link().\
                    using("g++").\
                    with_options([PYTHON_EXTENSION_LINKING_OPTIONS,OPENMP_OPTION]).\
                    using_libs([PYTHON_LIBRARY]).\
                    this_object_files([files_to_link["matrix"]]).\
                    to_produce("condensedMatrix.so")
    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')            
    os.system( linkDSL.getLinkingCommand())
    
    linkDSL = Link().\
                    using("g++").\
                    with_options([PYTHON_EXTENSION_LINKING_OPTIONS,OPENMP_OPTION]).\
                    using_libs([PYTHON_LIBRARY,"lapack","cudart"]).\
                    using_lib_locations([CUDA_LIB_64]).\
                    this_object_files([files_to_link["ThRMSDSerial"],files_to_link["ThRMSDSerialOmp"],files_to_link["ThRMSDCuda"],files_to_link["kernel_functions_serial"],\
                                       files_to_link["kernel_functions_omp"],files_to_link["kernel_functions_cuda"],files_to_link["NumpyHelperFuncs"],files_to_link["pyRMSD"],\
                                       files_to_link["RMSDomp"],files_to_link["RMSD"],files_to_link["RMSDSerial"],files_to_link["RMSDTools"],]).\
                    to_produce("pyRMSD_cfuncs.so")
    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
    os.system(linkDSL.getLinkingCommand())
    
    if options.use_cuda:
        linkDSL = Link().\
                        using("g++").\
                        with_options([]).\
                        using_libs(["cudart"]).\
                        using_lib_locations([CUDA_LIB_64]).\
                        this_object_files([files_to_link["ctest"],files_to_link["ThRMSDSerial"],files_to_link["kernel_functions_serial"],files_to_link["ThRMSDCuda"],\
                                           files_to_link["kernel_functions_cuda"],files_to_link["RMSD"]]).\
                        to_produce("test_main")
                        
        os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
        os.system(linkDSL.getLinkingCommand())
    
    os.system("mv pyRMSD_cfuncs.so pyRMSD/")
    os.system("mv condensedMatrix.so pyRMSD/")
    if options.use_cuda:
        os.system("mv test_main src/theobald/test/")
    
    ##Calculators
    if options.use_cuda:
        calcs_str = """
def availableCalculators():
    return {"PYTHON_CALCULATOR":-1,"SERIAL_CALCULATOR":0,"OMP_CALCULATOR":1,"THEOBALD_CUDA_CALCULATOR":2,"THEOBALD_SERIAL_CALCULATOR":3,"THEOBALD_SERIAL_OMP_CALCULATOR":4}
"""
    else:
        calcs_str = """
def availableCalculators():
    return {"PYTHON_CALCULATOR":-1,"SERIAL_CALCULATOR":0,"OMP_CALCULATOR":1,"THEOBALD_SERIAL_CALCULATOR":3,"THEOBALD_SERIAL_OMP_CALCULATOR":4}
"""
    os.system('echo "\033[33mWriting calculators...\033[0m"')
    open("pyRMSD/Calculators.py","w").write(calcs_str)
    
    os.system('echo "\033[32mCleaning...\033[0m"')
    for produced_file in  files_to_link:
        if files_to_link[produced_file] != "":
            os.system("rm "+files_to_link[produced_file])

