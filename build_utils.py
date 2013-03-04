import os

## Some really simple DSLs to make our life simpler (ok... just for fun!)
class Compile:
    def __init__(self):
        self.compiler =  ""
        self.compiler_options = ""
        self.includes = []
        self.file_name = ""
        
    def using(self, compiler):
        self.compiler = compiler
        return self
    
    def the_file(self, file_name):
        self.file_name = file_name
        return self
    
    def with_options(self, options):
        self.compiler_options = options
        return self
        
    def including_folder(self, folder):
        self.includes.append(folder)
        return self
    
    def including_folders(self, folders):
        self.includes.extend(folders)
        return self
    
    def getCompilingCommand(self):
        compiler_str = self.compiler + " "
        compiler_str += self.compiler_options + " "
        for include in self.includes:
            compiler_str+= "-I"+include + " "
        compiler_str+= "-c "+self.file_name
        return compiler_str

class Link:
    def __init__(self):
        self.linker =  ""
        self.options = []
        self.libs = []
        self.object_files = []
        self.product = ""
        self.lib_locations = []
        
    def using(self,linker):
        self.linker = linker
        return self
    
    def with_options(self,options):
        self.options.extend(options)
        return self
    
    def using_libs(self,libs):
        self.libs = libs
        return self
    
    def using_lib_locations(self,locations):
        self.lib_locations = locations
        return self
    
    def this_object_files(self, object_files):
        self.object_files = object_files
        return self
    
    def to_produce(self, product):
        self.product = product
        return self
    
    def getLinkingCommand(self):
        linking_str = self.linker + " "
        
        for option in self.options:
            linking_str += option + " "
        
        for o in self.object_files:
            linking_str += o+" "
        
        for lib_loc in self.lib_locations:
            linking_str += "-L"+lib_loc+" "
        
        for lib in self.libs:
            linking_str += "-l"+lib+" "
        
        linking_str += "-o "+self.product
        
        return linking_str

# And a very helpful function
def compile_a_file_collection(base_dir, file_collection, compiler, options, includes, product_extension, files_to_link):
    for folder in file_collection:
        os.chdir(folder)
        for file_name in file_collection[folder]:
            filewoext,extension = file_name.split(".") #@UnusedVariable
            files_to_link[filewoext] = folder+"/"+filewoext+product_extension
            comp = Compile().using(compiler).with_options(options).including_folders(includes).the_file(file_name)
            os.system('echo "\\033[31m'+ comp.getCompilingCommand()+'\\033[0m"')
            os.system( comp.getCompilingCommand())
        os.chdir(base_dir)

def get_object_file(all_objects, file_id):
    if not file_id in all_objects:
        return ""
    else:
        return all_objects[file_id]
    
