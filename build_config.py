import distutils.sysconfig
import numpy
import json
import os.path

# CUDA_BASE                # Base dir for CUDA installation
# CUDA_INCLUDE_FOLDER      # CUDA headers path
# CUDA_LIBRARIES_FOLDER    # CUDA libs path ( /lib if you're running it in a 32b machine)
# CUDA_ARCHITECHTURE       # CUDA architecture of your card.

def get_config_options_for(default_conf_file, this_conf_file):
    # Load defaults
    file_handler = open(default_conf_file,"r")
    base_conf = json.loads("".join(file_handler.readlines()))
    file_handler.close()
    
    # Override 
    if this_conf_file != None:
        file_handler = open(default_conf_file,"r")
        extra_conf = json.loads("".join(file_handler.readlines()))
        file_handler.close()
        for key in extra_conf:
            base_conf[key] = extra_conf["key"]
    
    # Process 
    if base_conf["PYTHON_INCLUDE_FOLDER"] == "AUTO":
        base_conf["PYTHON_INCLUDE_FOLDER"] = distutils.sysconfig.get_python_inc()
        
    if base_conf["PYTHON_LIBRARY_FOLDER"] == "AUTO":
        base_conf["PYTHON_LIBRARY_FOLDER"] = distutils.sysconfig.get_python_lib()
        
    if base_conf["NUMPY_INCLUDE"] == "AUTO":
        base_conf["NUMPY_INCLUDE"] = numpy.get_include()
        
    base_conf["CUDA_INCLUDE"] = os.path.join(base_conf["CUDA_BASE"],base_conf["CUDA_INCLUDE_FOLDER"])  
    base_conf["CUDA_LIBRARIES"] = os.path.join(base_conf["CUDA_BASE"],base_conf["CUDA_LIBRARIES_FOLDER"])  
    base_conf["CUDA_OPTIONS"] = base_conf["CUDA_OPTIONS"]%(base_conf["CUDA_ARCHITECHTURE"])
    
    base_conf["BASE_DIR"] = os.getcwd()
    
    return base_conf
    