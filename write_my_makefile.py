#!/usr/bin/python3
import argparse
from os import path, getcwd

parser = argparse.ArgumentParser()
parser.add_argument('directory', help="directory to write the makefile to")
parser.add_argument('program_name',help="name of the executable to generate (no path)")
parser.add_argument('--files','-f',help="C++ source files to compile.", nargs='+', type=str)
parser.add_argument('--PETSC_DIR','-PETSC_DIR', help="PETSc directory")
parser.add_argument('--PETSC_ARCH','-PETSC_ARCH',help="PETSc architecture")
parser.add_argument('--CXX','-CXX')
parser.add_argument('--CXXFLAGS','-CXXFLAGS')

ODEINT_DIR=str(path.abspath(getcwd()))
ODEINT_INCLUDE_DIR=ODEINT_DIR+'/include/'

def parse_args():
    args = {}
    pinfo = parser.parse_args()
    if pinfo.directory:
        args['outdir'] = path.abspath(pinfo.directory)
    if pinfo.program_name:
        args['target'] = pinfo.program_name
    if pinfo.files:
        args['srcs'] = [str(fl).rsplit('/',1)[1] for fl in pinfo.files]

        
    if pinfo.PETSC_DIR:
        args['PETSC_DIR'] = pinfo.PETSC_DIR
    else:
        args['PETSC_DIR'] = '/usr/local/petsc'
    if pinfo.PETSC_ARCH:
        args['PETSC_ARCH'] = pinfo.PETSC_ARCH
    else:
        args['PETSC_ARCH'] = 'arch-linux2-c-debug'
    if pinfo.CXX:
        args['CXX'] = pinfo.CXX
    else:
        args['CXX'] = 'clang++'
    if pinfo.CXXFLAGS:
        args['CXXFLAGS'] = pinfo.CXXFLAGS
    else:
        args['CXXFLAGS'] = '-std=c++17 -O3 -march=native -mtune=native -fopenmp'
    


    args['petsc_var_include'] = 'include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscvariables'
    args['petsc_rules_include'] = 'include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/petscrules'
    args['petsc_base_include'] = 'include $(PETSC_DIR)/$(PETSC_ARCH)/lib/petsc/conf/base'

    return args


def form_lines(args):
    lines = [f"CXX={args['CXX']}\n",
             f"CXXFLAGS={args['CXXFLAGS']}\n\n",
             f"PETSC_DIR={args['PETSC_DIR']}\n",
             f"PETSC_ARCH={args['PETSC_ARCH']}\n",
             f"{args['petsc_var_include']}\n",
             f"{args['petsc_rules_include']}\n",
             f"ODEINT_INCL={ODEINT_INCLUDE_DIR}\n",
             f"INCLS=-I$(ODEINT_INCL) $(PETSC_CC_INCLUDES) `pkg-config --cflags eigen3`\n\n",
             f"TARGET={args['target']}\n\n",
             f".PHONY: all clean $(TARGET)\n\n",
             f"all: $(TARGET)\n\n",
             f"$(TARGET): {' '.join(args['srcs'])}\n",
             f"\t$(CXX) $(CXXFLAGS) $^ -o $@ $(INCLS) $(PETSC_WITH_EXTERNAL_LIB)\n\n",
             f"clean:\n",
             f"\t@$(RM) *.o"]

    return lines
    

def write_makefile(args, lines):
    mfname = args['outdir']+'/makefile'
    with open(mfname, 'w') as makefile:
        for line in lines:
            makefile.write(line)
            




if __name__ == '__main__':
    args = parse_args()
    write_makefile(args, form_lines(args))
        

