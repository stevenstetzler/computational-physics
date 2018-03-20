import subprocess
import os
import numpy as np
from glob import glob

def slurm():
    header_lines = ['#!/bin/bash']

    out_file = '#SBATCH --output=metro-{0:0.1f}-{1:0.1f}.out'

    job_name = '#SBATCH --job-name="m {0:0.1f}-{1:0.1f}"'

    script_file = 'metro-{0:0.1f}-{1:0.1f}.sh'

    run_command = './metropolis {0} {1} {2} {3}'

    filename = 'metro-{0:0.1f}-{1:0.1f}.txt'

    num_scripts = 5
    #num_T = 125
    num_T = 12
    #Ts = np.linspace(0.01, 5, num_T)

    Ts = np.linspace(0.01, 0.452661, num_T)

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
        script_lines.append(out_file.format(T_low, T_high))

        script_lines.append(job_name.format(T_low, T_high))

        if high_idx != -1:
            script_lines.append(run_command.format(T_low, T_high, len(Ts[low_idx:high_idx + 1]), filename.format(T_low, T_high)))
        else:
            script_lines.append(run_command.format(T_low, T_high, len(Ts[low_idx:]), filename.format(T_low, T_high)))

        scripts.append(script_file.format(T_low, T_high))
        out = open(scripts[i], 'w')
        for l in script_lines:
            #print(l)
            out.write(l + '\n')
        out.close()
    
    #print(scripts)
    procs = []
    for script in scripts:
        procs.append(subprocess.Popen(['sbatch', script]))

def computer():
    code_files = glob("crit*")
    exec_files = [c for c in code_files if '.cpp' not in c]
    print(exec_files)

    filename = "crit_L_{0}_T_{1:0.1f}_{2:0.1f}.txt"
    command = "{0} {1} {2} {3} {4}"

    #num_T = 125
    num_T = 12

    num_runs = 1

    Ts = np.linspace(2, 3, num_T)

    for exec_file in exec_files:
        L = exec_file.split("_")[-1]

        for i in range(num_runs):

            low_idx = i * num_T / num_runs
            if i != 0:
                low_idx += 1
            
            T_low = Ts[low_idx]

            high_idx = (i + 1) * num_T / num_runs
            if i == num_runs - 1:
                high_idx = -1

            T_high = Ts[high_idx]
            
            #print(T_low, T_high)
            
            options = ['srun']
            options.append(exec_file)
            options.append(str(T_low))
            options.append(str(T_high))

            if high_idx != -1:
                options.append(str(len(Ts[low_idx:high_idx + 1])))
                # run_command = command.format(exec_file, T_low, T_high, len(Ts[low_idx:high_idx + 1]), filename.format(L, T_low, T_high))
            else:
                options.append(str(len(Ts[low_idx:])))
                # run_command = command.format(exec_file, T_low, T_high, len(Ts[low_idx:]), filename.format(L, T_low, T_high))
            
            options.append(filename.format(L, T_low, T_high))
            #print(options)

            options.append(">>")

            options.append("crit_L_{0}_T_{1}_{2}.out".format(L, T_low, T_high))

            options.append("&")

            run_command = " ".join(options)

            print(run_command)
            subprocess.Popen(options)

def main():
    computer()

if __name__ == "__main__":
    main()
