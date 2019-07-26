#!/bin/bash
SBATCH --account=def-jcraig1
SBATCH --mem-per-cpu=1024M
SBATCH --mail-user=erfan.amiri@uwaterloo.ca
SBATCH --mail-type=FAIL
SBATCH --time=02-00:00
SBATCH --job-name=SIMPER_REALIZATION_2

./SIMPER.exe 2