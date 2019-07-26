#!/bin/bash
# submit with:
#       sbatch submit-job-to-graham.sh     
SBATCH --account=def-jcraig1                      # your group NOT SURE IF THAT IS CALLED LIKE THIS
SBATCH --mem-per-cpu=1024M                        # memory; default unit is megabytes
SBATCH --mail-user=erfan.amiri@uwaterloo.ca       # email address for notifications
SBATCH --mail-type=FAIL                           # email send only in case of failure
SBATCH --time=0-00:15                             # time (DD-HH:MM);  here 15min
SBATCH --job-name=erfan-test-job                  # name of job in queque
# here all commands it needs to run your executable and do other stuff
# = things you would do one after each other on the comman line
./SIMPER.exe 12                  # job
