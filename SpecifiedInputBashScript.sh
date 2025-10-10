#!/bin/bash

# This bash script accompanies the TestPaperCommandLineVertexSimulationTest. This script will execute the test TestPaperCommandLineVertexSimulation.
# This script is set up to take in 10 specific varaibles for Lambda and Gamma from a CSV file.
# The script expects the output csv to only contain the values for Lambda and Gamma and the headers above each value :
# Note: this script expects that values are given in a comman delimited format
# Note this bash script should read all numbers in the bash script regardless of length
# So one could begin with a single set of parameters (low fidelity) then supply 10 in the following runs (high fidelity)

# ----------------   Begining of CSV -----------------
#       Lambda     Gamma      Runs    Simulation
#          1         10        3           1
#          2         9         2           2
#          3         8         2           3
#          4         7         1           4
#          5         6         1           5
#          6         5         3           6
#          7         4         2           7
#          8         3         1           8
#          9         2         1           9
#          10        1         3           10
# ------------------   End of CSV  --------------------

# Here we will create a loop that goes through each row and takes the given values for Lambda and Gamma

while IFS="," read -r rec_column1 rec_column2 rec_column3 rec_column4
do
i=1
## This for loop re-runs each parameter pair based on the Run number in the CSV
for ((i; i <= $rec_column3; i += 1)); do
     echo "Lambda: $rec_column1"
     echo "Gamma: $rec_column2"
     echo "Simulation: $rec_column4"
     echo "Run: $i"
     ## Running simulation with read in parameters and taking a random interger number between 1 - 1000 as the 5th command line option
     ## To prevent issues with the possibility of hanging simulations we introduce a watcher to check over the process and ensure it does not exceed a timeout.
     ## Currentlly this timeout is defaulted to 30 minutes (1800 seconds), which is about ten mintues over the expected runtime just incase a simulation is running slowly.
     ## If a process exceeds this limit it will be killed and a corresponding warning will be output
     ((~/build/projects/BayesianTissueProject2/test/TestPaperCommandLineSpeedSimulation --timeout 1 -opt1 $rec_column1 -opt2 $rec_column2 -opt3 $i -opt4 $rec_column4 -opt5 $((1+ $RANDOM % 1000))) & pid=$!
     (sleep 1800 && kill -HUP $pid) 2>/dev/null & watcher=$! 
     if wait $pid 2>/dev/null; then
        ##echo "The simulation ran correctly"
        pkill -HUP -P $watcher
        wait $watcher
     else 
        echo "The simulation was killed due to exceeding the timout limit"
        echo "The Parameter pair of Lambda = $rec_column1 and Gamma = $rec_column2 in Run $i failed due to timing out"
     fi)&      
     echo ""
done
done < <(tail -n  +1 ~/Chaste/projects/BayesianTissueProject2/ExampleCommandLineCSV.csv)
# Echos that the simulations have all been ran, this does not mean they did not error.
echo All done
