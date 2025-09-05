#!/bin/bash

# This bash script accompanies the TestPaperCommandLineVertexSimulationTest. This script will execute the test TestPaperCommandLineVertexSimulation.
# This script is set up to take in 10 specific varaibles for Lambda and Gamma from a CSV file.
# The script expects the output csv to only contain the values for Lambda and Gamma and the headers above each value :
# Note: this script expects that values are given in a comman delimited format
# Note this bash script should read all numbers in the bash script regardless of length
# So one could begin with a single set of parameters (low fidelity) then supply 10 in the following runs (high fidelity)

# ------------   Begining of CSV -------------
#       Lambda     Gamma      Runs
#          1         10        3
#          2         9         2
#          3         8         2
#          4         7         1
#          5         6         1
#          6         5         3
#          7         4         2
#          8         3         1
#          9         2         1
#          10        1         3
# --------------   End of CSV  ----------------

# Here we will create a loop that goes through each row and takes the given values for Lambda and Gamma


while IFS="," read -r rec_column1 rec_column2 rec_column3 rec_column4
do
## This for loop re-runs each parameter pair based on the Run number in the CSV
for ((i = 1; i <= $rec_column3; i += 1)); do
     echo "Lambda: $rec_column1"
     echo "Gamma: $rec_column2"
     echo "Simulation: $rec_column4"
     echo "Run: $i"
     ## Running simulation with read in parameters
     ~/bui/projects/BayesianTissueProject2/test/TestPaperCommandLineSpeedSimulation -opt1 $rec_column1 -opt2 $rec_column2 -opt3 $i -opt4 $rec_column4 &
     echo ""
done
done < <(tail -n +1 ~/Chaste/projects/BayesianTissueProject/ExampleCommandLineCSV.csv)
# Echos that the simulations have all been ran, this does not mean they did not error.
echo All done
