#!/bin/bash
#SBATCH --account=bgmp                    
#SBATCH --partition=bgmp   
#SBATCH --job-name=deduper
#SBATCH --output=%j.out
#SBATCH --error=%j.err                             
#SBATCH --cpus-per-task=8                 
#SBATCH --mem=32GB 

conda activate base

input_file=/projects/bgmp/tamho/bioinfo/Bi624/Deduper-tamdbho/input.sam
UMI_file=/projects/bgmp/tamho/bioinfo/Bi624/Deduper-tamdbho/STL96.txt
outputDIR=/projects/bgmp/tamho/bioinfo/Bi624/Deduper-tamdbho/output.test

/usr/bin/time -v ./ho_deduper.py -f $input_file -u $UMI_file -o $outputDIR
