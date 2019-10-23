from __future__ import print_function
from build_utils import compile_a_file_collection, Link
import optparse
import collections
import os
from build_config import get_config_options_for

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='%prog [--build-conf] [--cuda] [--build] [--clean] [--clean-all]', version='4.0.0')
    parser.add_option('--build-conf', dest = "conf_file",  help="Determines the file storing build configuration info.")
    parser.add_option('--cuda', dest = "cuda_type",  help="Use this flag if you want to compile the CUDA calculator.")
    parser.add_option('--build', dest = "build", action="store_true", help="Use this flag if you want to compile pyRMSD.")
    parser.add_option('--clean', dest = "clear",  action="store_true", help="Clear all .o generated files.")
    parser.add_option('--clean-all', dest = "clear_all", action="store_true",  help="The same as --clear, but it also removes generates libs and exes.")
    options, args = parser.parse_args()

    # Load configuration
    conf = get_config_options_for(os.path.join("build_conf","default.conf"), options.conf_file)


    ######################################
    ### SET CUDA FLAGS
    ######################################
    CUDA_PRECISION_FLAG = ""
    if not options.cuda_type is None:
        if options.cuda_type == "single":
            CUDA_PRECISION_FLAG = "CUDA_PRECISION_SINGLE"
        elif options.cuda_type =="double":
            CUDA_PRECISION_FLAG = "CUDA_PRECISION_DOUBLE"
        else:
            parser.error("Please, choose precision for CUDA building: 'single' or 'double'")
        options.use_cuda = True
        conf["CUDA_OPTIONS"] = conf["CUDA_OPTIONS"] +" -D"+CUDA_PRECISION_FLAG
        conf["DEFINE_USE_CUDA"] = conf["DEFINE_USE_CUDA"] +" -D"+CUDA_PRECISION_FLAG
    else:
        options.use_cuda = False
    ######################################

    ######################################
    ### FILE DESCRIPTION
    ######################################
    files_to_compile_with_nvcc = {
                                  "src/calculators/QCP":["QCPCUDAKernel.cu","QCPCUDAMemKernel.cu","kernel_functions_cuda.cu"]
    }

    files_to_compile_with_gcc = {
                                 "src/calculators":["RMSDCalculator.cpp","RMSDTools.cpp","RMSDCalculationData.cpp","KernelFunctions.cpp"],
                                 "src/calculators/factory":["RMSDCalculatorFactory.cpp"],
                                 "src/calculators/KABSCH":["KABSCHSerialKernel.cpp"],
                                 "src/calculators/QTRFIT":["QTRFITSerialKernel.cpp"],
                                 "src/calculators/QCP":["QCPSerialKernel.cpp","QCPSerialFloatKernel.cpp"],
                                 "src/calculators/NOSUP":["NOSUPSerialKernel.cpp"],
                                 "src/matrix":["Matrix.cpp","Statistics.cpp"],
                                 "src/python":["pyRMSD.cpp"],
                                 "src/pdbreaderlite":["PDBReader.cpp","PDBReaderObject.cpp"],
                                 "src/calculators/test":["main.cpp","test_tools.cpp","tests.cpp"],
    }

    files_to_compile_with_gcc_and_openmp = {
                                            "src/calculators/KABSCH":["KABSCHOmpKernel.cpp"],
                                            "src/calculators/QTRFIT":["QTRFITOmpKernel.cpp"],
                                            "src/calculators/QCP":["QCPOmpKernel.cpp"],
                                            "src/calculators/NOSUP":["NOSUPOmpKernel.cpp"]
    }
    ######################################

    #########################################
    ###
    ###     BUILDING PROCESS
    ###
    #########################################
    if options.build:
        files_to_link = collections.defaultdict(str)

        if options.use_cuda:
            conf["PYTHON_EXTENSION_OPTIONS"] = conf["PYTHON_EXTENSION_OPTIONS"]+" "+conf["DEFINE_USE_CUDA"]

        compile_a_file_collection(conf["BASE_DIR"], files_to_compile_with_gcc, "gcc", conf["PYTHON_EXTENSION_OPTIONS"], [conf["PYTHON_INCLUDE_FOLDER"], conf["NUMPY_INCLUDE"]], ".o",files_to_link)

        compile_a_file_collection(conf["BASE_DIR"], files_to_compile_with_gcc_and_openmp, "gcc", conf["PYTHON_EXTENSION_OPTIONS"]+" "+conf["OPENMP_OPTION"], [conf["PYTHON_INCLUDE_FOLDER"], conf["NUMPY_INCLUDE"]], ".o",files_to_link)

        if options.use_cuda:
            compile_a_file_collection(conf["BASE_DIR"], files_to_compile_with_nvcc, "nvcc", conf["CUDA_OPTIONS"], [conf["CUDA_INCLUDE"]], ".o",files_to_link)

        linkDSL = Link().\
                        using("g++").\
                        with_options([conf["PYTHON_EXTENSION_LINKING_OPTIONS"],conf["OPENMP_OPTION"]]).\
                        using_libs([conf["PYTHON_LIBRARY"]]).\
                        using_lib_locations([conf["PYTHON_LIBRARY_FOLDER"]]).\
                        this_object_files([files_to_link["Matrix"],files_to_link["Statistics"]]).\
                        to_produce("condensedMatrix.so")

        os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
        os.system( linkDSL.getLinkingCommand())

        linkDSL = Link().\
                        using("g++").\
                        with_options([conf["PYTHON_EXTENSION_LINKING_OPTIONS"],conf["OPENMP_OPTION"]]).\
                        using_libs([conf["PYTHON_LIBRARY"]]).\
                        using_lib_locations([conf["PYTHON_LIBRARY_FOLDER"]]).\
                        this_object_files([files_to_link["PDBReaderObject"],files_to_link["PDBReader"]]).\
                        to_produce("pdbReader.so")

        os.system('echo "\033[34m'+ linkDSL.getLinkingCommand()+'\033[0m"')
        os.system( linkDSL.getLinkingCommand())

        calculator_obj_files = [
                                files_to_link["RMSDTools"],
                                files_to_link["KernelFunctions"],
                                files_to_link["KABSCHSerialKernel"],
                                files_to_link["KABSCHOmpKernel"],
                                files_to_link["QTRFITSerialKernel"],
                                files_to_link["QTRFITOmpKernel"],
                                files_to_link["QCPSerialKernel"],
                                files_to_link["QCPSerialFloatKernel"],
                                files_to_link["QCPOmpKernel"],
                                files_to_link["NOSUPSerialKernel"],
                                files_to_link["NOSUPOmpKernel"],
                                files_to_link["RMSDCalculatorFactory"],
                                files_to_link["RMSDCalculationData"],
                                files_to_link["RMSDCalculator"],
                                files_to_link["pyRMSD"]
        ]

        if options.use_cuda:
            calculator_obj_files.extend([files_to_link["QCPCUDAKernel"],files_to_link["QCPCUDAMemKernel"],files_to_link["kernel_functions_cuda"]])
            calculator_libraries = [conf["PYTHON_LIBRARY"],conf["CUDA_LIBRARY"]]
            calculator_library_locations  = [conf["PYTHON_LIBRARY_FOLDER"], conf["CUDA_LIBRARIES"]]

        else:
            calculator_libraries = [conf["PYTHON_LIBRARY"]]
            calculator_library_locations  = [conf["PYTHON_LIBRARY_FOLDER"]]

        linkDSL = Link().\
                        using("g++").\
                        with_options([conf["PYTHON_EXTENSION_LINKING_OPTIONS"],conf["OPENMP_OPTION"]]).\
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
                        with_options([conf["OPENMP_OPTION"]]).\
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
            #"QCP_SERIAL_FLOAT_CALCULATOR":7,
            "QCP_OMP_CALCULATOR":8,
            "QCP_CUDA_CALCULATOR":9,
            "QCP_CUDA_MEM_CALCULATOR":10,
            "NOSUP_SERIAL_CALCULATOR":11,
            "NOSUP_OMP_CALCULATOR":12,
            #"NOSUP_CUDA_CALCULATOR":13
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
            #"QCP_SERIAL_FLOAT_CALCULATOR":7,
            "QCP_OMP_CALCULATOR":8,
            #"QCP_CUDA_CALCULATOR":9,
            #"QCP_CUDA_MEM_CALCULATOR":10,
            "NOSUP_SERIAL_CALCULATOR":11,
            "NOSUP_OMP_CALCULATOR":12,
            #"NOSUP_CUDA_CALCULATOR":13
    }
"""
        os.system('echo "\033[33mWriting available calculators...\033[0m"')
        open("pyRMSD/availableCalculators.py","w").write(calcs_str)

        # Save all produced files
        produced_file_handler = open(".products","w")
        for produced_file in  files_to_link:
            if files_to_link[produced_file] != "":
                produced_file_handler.write(files_to_link[produced_file] +"\n")
        produced_file_handler.close()
    ######################################

    #########################################
    ###
    ###     REMOVE ALL .o AND TRACKING FILES
    ###
    #########################################
    if options.clear or options.clear_all:
        os.system('echo "\033[32mCleaning...\033[0m"')
        if os.path.exists(".products"):
            produced_file = open(".products","r")
            for produced_file_line in produced_file:
                os.system("rm "+produced_file_line)
                print("rm "+produced_file_line[:-1])
            produced_file.close()
            # Clear the products file itself
            os.system("rm .products")
        # Remove all trackers
        os.system("find src/ -name '.modif*' -exec rm {} \;")
        # remove .pyc
        os.system("find src/ -name '*.pyc' -exec rm {} \;")
    ######################################

    #########################################
    ###
    ###     REMOVE ALL LIBS AND EXES
    ###
    #########################################
    if options.clear_all:
        os.system('echo "\033[32mCleaning exes and libs...\033[0m"')
        os.system("rm pyRMSD/calculators.so")
        os.system("rm pyRMSD/condensedMatrix.so")
        os.system("rm pyRMSD/pdbReader.so")
        os.system("rm src/calculators/test/test_rmsdtools_main")
    ######################################
