import subprocess
import os
import numpy as np
from glob import glob

def slurm():
    for L in [30, 60, 120, 240, 360]:
        header_lines = ['#!/bin/bash']

        out_file = '#SBATCH --output=crit-L-{0}-{1:0.2f}-{2:0.2f}.out'

        job_name = '#SBATCH --job-name="cL{0}-{1:0.2f}-{2:0.2f}"'

        script_file = 'crit-L-{0}-{1:0.2f}-{2:0.2f}.sh'

        run_command = './crit_L_{0} {1} {2} {3} {4}'

        filename = 'crit-L-{0}-{1:0.2f}-{2:0.2f}.txt'

        num_scripts = 15
        #num_T = 125
        num_T = 200
        Ts = np.linspace(1.26918531, 3.26918531, num_T)

        #Ts = np.linspace(0.01, 0.452661, num_T)

        scripts = []
        #print(Ts)
        for i in range(num_scripts):

            low_idx = i * num_T / num_scripts
            if i != 0:
                low_idx += 1
            
            T_low = Ts[low_idx]

            high_idx = (i + 1) * num_T / num_scripts
            if i == num_scripts - 1:
                high_idx = -1

            T_high = Ts[high_idx]
            
            #print(T_low, T_high)

            script_lines = [l for l in header_lines]
            script_lines.append(out_file.format(L, T_low, T_high))

            script_lines.append(job_name.format(L, T_low, T_high))

            if high_idx != -1:
                script_lines.append(run_command.format(L, T_low, T_high, len(Ts[low_idx:high_idx + 1]), filename.format(L, T_low, T_high)))
            else:
                script_lines.append(run_command.format(L, T_low, T_high, len(Ts[low_idx:]), filename.format(L, T_low, T_high)))

            scripts.append(script_file.format(L, T_low, T_high))
            out = open(scripts[i], 'w')
            for l in script_lines:
                print(l)
                out.write(l + '\n')
            out.close()
        
        #print(scripts)
        procs = []
        for script in scripts:
            procs.append(subprocess.Popen(['sbatch', script]))

def compile(scripts):
    procs = []

    for script in scripts:
        command = ['g++', '-O1', '-O2', '-O3']
        
        exec_name = script.split(".")[0]
        
        command.append('-o{0}'.format(exec_name))
        command.append(script)

        procs.append(subprocess.Popen(command))

    print(procs)

def computer():
    code_files = glob("crit*.cpp")
    compile(code_files)


def main():
    slurm()

if __name__ == "__main__":
    main()
