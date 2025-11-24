#!/bin/bash

# This bash script accompanies the TestPaperCommandLineVertexSimulationTest. This script will execute the test TestPaperCommandLineVertexSimulation.
# This script is set up to take in 10 specific varaibles for Lambda and Gamma from a CSV file.
# The script expects the output csv to only contain the values for Lambda and Gamma and the headers above each value :
# Note: this script expects that values are given in a comman delimited format
# Note this bash script should read all numbers in the bash script regardless of length
# So one could begin with a single set of parameters (low fidelity) then supply 10 in the following runs (high fidelity)

# ---------------------   Begining of CSV ----------------------
#       Lambda     Gamma      Runs    Simulation    TimeStep  
#          1         10        3           1          0.01
#          2         9         2           2          0.01
#          3         8         2           3          0.01
#          4         7         1           4          0.01
#          5         6         1           5          0.01
#          6         5         3           6          0.01
#          7         4         2           7          0.01
#          8         3         1           8          0.01
#          9         2         1           9          0.01
#          10        1         3           10         0.01
# -----------------------   End of CSV  -------------------------

# Here we will create a loop that goes through each row and takes the given values for Lambda and Gamma

while IFS="," read -r rec_column1 rec_column2 rec_column3 rec_column4 rec_column5
do
i=1
## This for loop re-runs each parameter pair based on the Run number in the CSV
for ((i; i <= $rec_column3; i += 1)); do
     echo "Lambda: $rec_column1"
     echo "Gamma: $rec_column2"
     echo "Simulation: $rec_column4"
     echo "Run: $i"
     echo "Time Step $rec_column5"
     ## Running simulation with read in parameters and taking a random interger number between 1 - 1000 as the 5th command line option
     ## To prevent issues with the possibility of hanging simulations we introduce a watcher to check over the processes output file .
     ## If the given outut file is not incrementign we assume the test is hanging and kill its pid
     ## Currentlly this timeout is defaulted to 5 minutes (300 seconds).
     ## If a process exceeds this limit it will be killed and a corresponding warning will be output
     ## WARNING: If you increase the time step above the default value of 0.01 you will need to increase the watchers tolerance

     ## If your test files are output somewhere other than /tmp/$USER/testoutput then you will need to change everything before /TestBayesianCommandLineRun2 to account for this
     
     (~/build/projects/BayesianTissueProject2/test/TestPaperCommandLineSpeedSimulation --timeout 1 -opt1 $rec_column1 -opt2 $rec_column2 -opt3 $i -opt4 $rec_column4 -opt5 $((1+ $RANDOM % 1000)) -opt6 $rec_column5) & pid=$!
      ./checkusage.sh "/tmp/"$USER"/testoutput/TestBayesianCommandLineRun2/_Sim_Number_"$rec_column4"_Lambda__"$rec_column1"_Gamma_"$rec_column2"_Run_"$i"" $pid

     echo ""
done
done < <(tail -n  +1 ~/Chaste/projects/BayesianTissueProject2/ExampleCommandLineCSV.csv)
# Echos that the simulations have all been ran, this does not mean they did not error.
echo All done