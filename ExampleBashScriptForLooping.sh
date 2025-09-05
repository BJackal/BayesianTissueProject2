#!/bin/bash

# This bash script accompanies the TestPaperCommandLineVertexSimulationTest. This script will execute the test  TestPaperCommandLineVertexSimulation.
# This example script shows how one could loop over varaibles to a given number and run simulations for each.
# Here we will declare some values we wish to later pass to a for loop.

counter=0
N=10
L=10
# The counter will echo how many times the N varaibles has been run
echo $counter
# Echo the current variable being run
echo "Current Run Lambda: $1 and Gamma: $2"
# This for loop will give us 10
for ((i = 0; i <= N; i += 1)); do
for ((j = 1; j <= L; j += 1)); do
	~/build/projects/BayesianTissueProject2/test/TestPaperCommandLineSpeedSimulation -opt1 $i -opt2 $j &
done
done
echo All done
