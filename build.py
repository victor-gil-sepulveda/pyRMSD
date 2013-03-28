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
PYTHON_EXTENSION_OPTIONS = "-pthread -g -fno-strict-aliasing -fmessage-length=0 -O3 -Wall \
-D_FORTIFY_SOURCE=2 -fstack-protector -funwind-tables -fasynchronous-unwind-tables -fPIC"
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
    parser.add_option('--cuda', dest = "cuda_type",  help="Use this flag if you want to compile the CUDA calculator.")
    options, args = parser.parse_args()
    
    CUDA_PRECISION_FLAG = ""
    if not options.cuda_type is None:
        if options.cuda_type == "single":
            CUDA_PRECISION_FLAG = "CUDA_PRECISION_SINGLE"
        elif options.cuda_type =="double":
            CUDA_PRECISION_FLAG = "CUDA_PRECISION_DOUBLE"
        else:
            parser.error("Please, choose  precision for CUDA building: 'single' or 'double'")
        options.use_cuda = True
        CUDA_OPTIONS = CUDA_OPTIONS +" -D"+CUDA_PRECISION_FLAG
        DEFINE_USE_CUDA = DEFINE_USE_CUDA +" -D"+CUDA_PRECISION_FLAG
    else:
        options.use_cuda = False
    
    ######################################
    ### Files we are going to use
    ######################################
    files_to_compile_with_nvcc = {
                                  "src/calculators/QCP":["QCPCUDAKernel.cu","kernel_functions_cuda.cu"]
    }
    
    files_to_compile_with_gcc = {
                                 "src/calculators":["RMSDCalculator.cpp","RMSDTools.cpp"],
                                 "src/calculators/factory":["RMSDCalculatorFactory.cpp"],
                                 "src/calculators/KABSCH":["KABSCHSerialKernel.cpp"],
                                 "src/calculators/QTRFIT":["QTRFITSerialKernel.cpp"],
                                 "src/calculators/QCP":["QCPSerialKernel.cpp"],
                                 "src/matrix":["Matrix.cpp","Statistics.cpp"],
                                 "src/python":["pyRMSD.cpp","readerLite.cpp"],
                                 "src/pdbreaderlite":["PDBReader.cpp"],
                                 "src/calculators/test":["main.cpp","test_tools.cpp","tests.cpp"],
    }
    
    files_to_compile_with_gcc_and_openmp = {
                                            "src/calculators/KABSCH":["KABSCHOmpKernel.cpp"],
                                            "src/calculators/QTRFIT":["QTRFITOmpKernel.cpp"],
                                            "src/calculators/QCP":["QCPOmpKernel.cpp"],
    }
    #########################################
    
    files_to_link = collections.defaultdict(str)
    
    if options.use_cuda:
        PYTHON_EXTENSION_OPTIONS = PYTHON_EXTENSION_OPTIONS+" "+DEFINE_USE_CUDA
    
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
    
    
    calculator_obj_files = [
                            files_to_link["RMSDTools"],
                            files_to_link["KABSCHSerialKernel"],
                            files_to_link["KABSCHOmpKernel"],
                            files_to_link["QTRFITSerialKernel"],
                            files_to_link["QTRFITOmpKernel"],
                            files_to_link["QCPSerialKernel"],
                            files_to_link["QCPOmpKernel"],
                            files_to_link["RMSDCalculatorFactory"],
                            files_to_link["RMSDCalculator"],
                            files_to_link["pyRMSD"]
    ]
    if options.use_cuda:
        calculator_obj_files.extend([files_to_link["QCPCUDAKernel"],files_to_link["kernel_functions_cuda"]])
        calculator_libraries = [PYTHON_LIBRARY,CUDA_LIBRARY]
        calculator_library_locations  = [PYTHON_LIBRARY_FOLDER, CUDA_LIBRARIES_FOLDER]
        
    else:
        calculator_libraries = [PYTHON_LIBRARY]
        calculator_library_locations  = [PYTHON_LIBRARY_FOLDER]
    
    linkDSL = Link().\
                    using("g++").\
                    with_options([PYTHON_EXTENSION_LINKING_OPTIONS,OPENMP_OPTION]).\
                    using_libs(calculator_libraries).\
                    using_lib_locations(calculator_library_locations).\
                    this_object_files(calculator_obj_files).\
                    to_produce("calculators.so")
    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
    os.system(linkDSL.getLinkingCommand())
    
    test_obj_files = list(calculator_obj_files)
    test_obj_files.remove(files_to_link["pyRMSD"])
    test_obj_files.extend([files_to_link["main"], files_to_link["test_tools"], files_to_link["tests"]])
    linkDSL = Link().\
                    using("g++").\
                    with_options([OPENMP_OPTION]).\
                    using_libs(calculator_libraries).\
                    using_lib_locations(calculator_library_locations).\
                    this_object_files(test_obj_files).\
                    to_produce("test_rmsdtools_main")
                    
    os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
    os.system(linkDSL.getLinkingCommand())
    
    os.system("mv calculators.so pyRMSD/")
    os.system("mv condensedMatrix.so pyRMSD/")
    os.system("mv pdbReader.so pyRMSD/")
    os.system("mv test_rmsdtools_main src/calculators/test")
    
    ##Calculators
    if options.use_cuda:
        calcs_str = """
def availableCalculators():
    return {
            "KABSCH_SERIAL_CALCULATOR": 0, 
            "KABSCH_OMP_CALCULATOR":1, 
            #"KABSCH_CUDA_CALCULATOR":2, 
            "QTRFIT_SERIAL_CALCULATOR":3,
            "QTRFIT_OMP_CALCULATOR":4,
            #"QTRFIT_CUDA_CALCULATOR":5,
            "QCP_SERIAL_CALCULATOR":6,
            "QCP_OMP_CALCULATOR":7,
            "QCP_CUDA_CALCULATOR":8
    }
"""
    else:
        calcs_str = """
def availableCalculators():
    return {
            "KABSCH_SERIAL_CALCULATOR": 0, 
            "KABSCH_OMP_CALCULATOR":1, 
            #"KABSCH_CUDA_CALCULATOR":2, 
            "QTRFIT_SERIAL_CALCULATOR":3,
            "QTRFIT_OMP_CALCULATOR":4,
            #"QTRFIT_CUDA_CALCULATOR":5,
            "QCP_SERIAL_CALCULATOR":6,
            "QCP_OMP_CALCULATOR":7,
            #"QCP_CUDA_CALCULATOR":8
    }
"""
    os.system('echo "\033[33mWriting available calculators...\033[0m"')
    open("pyRMSD/availableCalculators.py","w").write(calcs_str)
    
    os.system('echo "\033[32mCleaning...\033[0m"')
    for produced_file in  files_to_link:
        if files_to_link[produced_file] != "":
            os.system("rm "+files_to_link[produced_file])

