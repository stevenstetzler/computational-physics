import subprocess
import os
import numpy as np
from glob import glob

def compile(scripts):
    procs = []

    for script in scripts:
        command = ['g++', '-O1', '-O2', '-O3']
        
        exec_name = script.split(".")[0]
        
        command.append('-o{0}'.format(exec_name))
        command.append(script)

        procs.append(subprocess.Popen(command))

    print(procs)

def main():
    code = glob("crit*.cpp")
    compile(code)

if __name__ == "__main__":
    main()
