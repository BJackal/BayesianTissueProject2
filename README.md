# A Chaste User Project for inputing specific parameters from CSV into a Chaste simulation at command line

This Chaste project contains source code necessary to run simulations and generate results in line with the work of Kursawe et al (https://doi.org/10.1016/j.jtbi.2018.01.020).
Generating multiple vertex model based tissues from the command line based on input perimeter contractility and line tension parameters.

This project contains a base test which can simulate tissues based on input hard coded parameters and a second test which takes parameters from a csv file and inputs them as command line arguments.

**Testing requirements**
- [x] Working Base Test
- [x] Functionality for test to interact with command line inputs
    
**Writer outputs from Chaste that are necessary**

- [x] Polygon Number
- [x] Area Ratios - No explicit writer. Can be performed through combination of the polygon and Area writers.
- [x] Cell perimeter
- [x] Edge length - Note that this only produces the edge lengths for internal edges and not those on the boundary
- [x] Cell elongation
- [x] Area deviation
- [x] Area correlation
- [x] Polygon number correlation
- [x] Neighbour Areas
- [x] Neighbour numbers
- [x] Neighbour numbers correlation

Need looking at further
- [ ] Laser recoil
- [ ] Area asymetry
- [ ] Perimeter asymetry

# Instructions for running the project
**Installing Chaste**
This project is designed and intended to be used as a Chaste user project.
So first you will need to download and install Chaste and its dependencies. A guide for doing this can be found here: https://chaste.github.io/docs/

**Dowloading the project**
After installing Chaste the simplest way to dowloading this project as a zip file and exporting it to the Chaste project folder (PATH_TO_CHASTE/Chaste/projects).
You will then need to follow the steps outlined at https://chaste.github.io/docs/user-guides/user-projects/ for installing,setting up and running the project.

**Change corresponding files in trunk**
To allow for some of the code within this project to run we are required to change VertexBasedCellPopulation and MutableVertexMesh in the main Chaste trunk code. This is to move some functions from private to protected such that child classes can access them.

To do the new versions of the files in ChangeCorresponding trunk and replace *both hpp and cpp* files of VertexBasedCellPopulation (trunk address: Chaste/cell_based/src/population ) and MutableVertexMesh (trunk address: Chaste/mesh/src/vertex).

**Running the tests**

This project contains two main tests. 

First, TestPaperVertexSimulation.hpp which is a standard version of a chaste test and can be run by typing "ctest -V -R TestPaperSpeedSimulation".
This will execute the test with whatever parameters are currently in the file.

Second, TestPaperCommandLineSpeedSimulation.hpp. This test cannot be run using ctest as it requires inputs from the command line. For ease a bash script can be found in BayesianTissueProject/ExampleBashScriptForLooping.sh and a example csv.
To run this test first follow the steps for setting up the the project (this assumes you are in your build directory)then run "ccmake PATH_TO_CHASTE/Chaste && make -j4 projects/BayesianTissueProject && cd  /PATH_TO_CHASTE/Chaste/projects/BayesianTissuePorjects && bash ./ExampleBashScriptForLooping.sh". This will execute the example test with the example csv file. 

To change the input csv file simply change the file name in the bash script and palce your new target csv in the same folder. Note, if you simply run the test with a new csv and no code is changed you simply need to change the bash script then run  
"bash ./ExampleBashScriptForLooping.sh".
