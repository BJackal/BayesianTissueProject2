/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTPAPERSPEEDSIMULATION_HPP_
#define TESTPAPERSPEEDSIMULATION_HPP_

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "HoneycombVertexMeshGenerator.hpp"
#include "VoronoiVertexMeshGenerator.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include "TargetAreaLinearGrowthModifier.hpp"
#include "FarhadifarForce.hpp"
#include "SimpleTargetAreaModifier.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "NoCellCycleModel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "FixedSequenceCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "TransitCellProliferativeType.hpp"
#include "CellVolumesWriter.hpp"
#include "CellAgesWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellDataItemWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellMutationStatesCountWriter.hpp"
#include "CellPropertyRegistry.hpp"
#include "SmartPointers.hpp"
#include "WildTypeCellMutationState.hpp"
#include "StemCellProliferativeType.hpp"

#include "CellForcesWriter.hpp"
#include "VertexModelDataWriter.hpp"
#include "CellEdgeCountWriter.hpp"
#include "CellPerimeterWriter.hpp"
#include "ForwardEulerNumericalMethod.hpp"

#include "AreaCorrelationWriter.hpp"
#include "NeighbourNumberCorrelationWriter.hpp"
#include "VertexEdgeLengthWriter.hpp"

#include "TargetAreaLinearGrowthModifier.hpp"
#include "CellForcesWriter.hpp"
#include "RK4NumericalMethod.hpp"
#include "TransitCellProliferativeType.hpp"
#include "ModifiedVertexBasedCellPopulation.hpp"
#include "CellsGenerator.hpp"
#include "FixedSequenceCellCycleModel.hpp"
#include "CellCycleTimesGenerator.hpp"
#include "ExtendedHoneycombVertexMeshGenerator.hpp"
#include "FasterMutableVertexMesh.hpp"
#include "VertexModelDataWriter.hpp"
#include "FarhadifarForceWriter.hpp"
#include "PolygonNumberCorrelationWriter.hpp"

#include "CommandLineArguments.hpp"


/**
 * These tests check and demonstrate simulation of vertex based models with edge  Srn models
 */
class TestPaperSpeedSimulation : public AbstractCellBasedTestSuite
{

public:

    /*
     * Test running vertex based model with edge based SRN.
     */
    void TestRunningBayesianTissueSpeed()
    {
        int mRandomSeed = 0;
        //double mDt = 0.01;
        unsigned mNumberGenerations = 7u; //o
        double mAverageCellCycleTime = 20.0; // o
        //double mEndTime = 4.0*7.0*100.0;
        double mNewEdgeLengthFactor = 1.5; // o
        bool mRestrictVertexMovement = true; // o
        bool mRandomiseT1SwapOrder = false;
        //unsigned mSamplingTimestepMultiple = 0u;
        double mT1SwapThreshold = 0.01;   // o
        double mT2SwapThreshold = 0.001;  // o
        unsigned mInitialSize = 2u; // o
        double mLineTensionParameter = 0.12; // o // Lambda -0.85, 0.0 , 0.12
        double mPerimeterContractilityParameter = 0.04; // o   // Gamma 0.1 , 0.1 , 0.04
        double mBoundaryTensionParameter = mLineTensionParameter; // o
        //bool mUseRungeKuttaMethod = false; //o

        CellCycleTimesGenerator* p_cell_cycle_times_generator = CellCycleTimesGenerator::Instance();

        if(mRandomSeed == 0u)
     	{
	      p_cell_cycle_times_generator->SetRandomSeed( time(NULL) );
	    }
	    else
	    {
		  p_cell_cycle_times_generator->SetRandomSeed( mRandomSeed );
    	}
        p_cell_cycle_times_generator->SetRate(3.0/(2.0*mAverageCellCycleTime));
        p_cell_cycle_times_generator->GenerateCellCycleTimeSequence();

        /* First we create a regular vertex mesh. */
        ExtendedHoneycombVertexMeshGenerator generator(mInitialSize, mInitialSize, false, mT1SwapThreshold, mT2SwapThreshold, 1.0);
        FasterMutableVertexMesh<2, 2>* p_mesh = generator.GetMesh();
        p_mesh->SetRandomizeT1SwapOrderBoolean( mRandomiseT1SwapOrder );

        p_mesh->SetCellRearrangementRatio(mNewEdgeLengthFactor);
        p_mesh->SetCheckForInternalIntersections(true);


        std::vector<CellPtr> cells;
        MAKE_PTR(WildTypeCellMutationState, p_state);
        MAKE_PTR(TransitCellProliferativeType, p_diff_type);

        for (unsigned elem_index=0; elem_index < p_mesh->GetNumElements(); elem_index++)
        {
            /* Initalise cell cycle */
            FixedSequenceCellCycleModel* p_cc_model = new FixedSequenceCellCycleModel();
            //UniformG1GenerationalCellCycleModel* p_cc_model = new UniformG1GenerationalCellCycleModel();
            p_cc_model->SetDimension(2);
            p_cc_model->SetG2Duration(1.0/3.0*mAverageCellCycleTime);
            p_cc_model->SetMDuration(1e-12);
            p_cc_model->SetSDuration(1e-12);
            p_cc_model->SetMaxTransitGenerations(mNumberGenerations);

            CellPtr p_cell(new Cell(p_state, p_cc_model));
            p_cell->SetCellProliferativeType(p_diff_type);


            //double birth_time = -RandomNumberGenerator::Instance()->ranf()*12.0;
            double birth_time = 0.0;
            p_cell->SetBirthTime(birth_time);
            cells.push_back(p_cell);
        }

        /* Using the vertex mesh and cells, we create a cell-based population object, and specify which results to
         * output to file. */
        ModifiedVertexBasedCellPopulation<2> cell_population(*p_mesh, cells);

        cell_population.SetRestrictVertexMovementBoolean(mRestrictVertexMovement);


        // Cell writers
        cell_population.AddCellWriter<VertexModelDataWriter>();
        //cell_population.AddCellWriter<CellProliferativePhasesWriter>();
        cell_population.AddCellWriter<CellAgesWriter>();
        cell_population.AddCellWriter<CellEdgeCountWriter>();
        cell_population.AddCellWriter<CellPerimeterWriter>();

        // Cell Population Writers
        cell_population.AddCellPopulationCountWriter<FarhadifarForceWriter>(); //TODO
        cell_population.AddCellPopulationCountWriter<AreaCorrelationWriter>();
        cell_population.AddCellPopulationCountWriter<PolygonNumberCorrelationWriter>(); //TODO
        cell_population.AddCellPopulationCountWriter<NeighbourNumberCorrelationWriter>();
        cell_population.AddPopulationWriter<VertexEdgeLengthWriter>(); 

        //cell_population.rGetMesh().SetCellRearrangementThreshold(0.2);

        cell_population.SetWriteCellVtkResults(false);
        cell_population.SetWriteEdgeVtkResults(false);

        CellPtr p_cell_0b = cell_population.GetCellUsingLocationIndex(3);

        double initial_target_area = 1.0;


        for (typename AbstractCellPopulation<2>::Iterator cell_iter = cell_population.Begin();
                cell_iter != cell_population.End();
                ++cell_iter)
        {
            // target areas
            cell_iter->GetCellData()->SetItem("target area", initial_target_area);
        }

        //double actual_end_time = (mAverageCellCycleTime * mNumberGenerations)*4.0;
        double actual_end_time = 700.0;

        /* We are now in a position to create and configure the cell-based simulation object, pass a force law to it,
         * and run the simulation. */
        OffLatticeSimulation<2> simulator(cell_population);
        simulator.SetOutputDirectory("TestBaysianSimulation");
        simulator.SetSamplingTimestepMultiple(100);
        simulator.SetDt(0.01);
        simulator.SetEndTime(actual_end_time);

        MAKE_PTR(FarhadifarForce<2>, p_force);

        p_force->SetPerimeterContractilityParameter(mPerimeterContractilityParameter); // Gamma 0.1 , 0.1 , 0.04
        p_force->SetLineTensionParameter(mLineTensionParameter); // Lambda -0.85, 0.0 , 0.12
        // If our Line tension parameter is negative in the bulk we need to set it to zero at the boundary
        // This is to prevent non-physical behaviour from occuring
        if(p_force->GetLineTensionParameter() < 0){
                    p_force->SetBoundaryLineTensionParameter(0.0);
        } else {
            p_force->SetBoundaryLineTensionParameter(mBoundaryTensionParameter);
        }
        //MAKE_PTR(NagaiHondaForce<2>, p_force);
        //MAKE_PTR(WelikyOsterForce<2>, p_force);
        simulator.AddForce(p_force);

        MAKE_PTR(TargetAreaLinearGrowthModifier<2>, p_growth_modifier);
        simulator.AddSimulationModifier(p_growth_modifier);

        // Pass an adaptive numerical method to the simulation
        boost::shared_ptr<AbstractNumericalMethod<2,2> > p_method(new ForwardEulerNumericalMethod<2,2>());
        p_method->SetUseAdaptiveTimestep(true);
        simulator.SetNumericalMethod(p_method);

        //MAKE_PTR(VertexBoundaryRefinementModifier<2>, refinement_modifier);
        //refinement_modifier->SetMaxEdgeLength(2.0);
        //refinement_modifier->SetMinEdgeLength(0.1);
        //simulator.AddSimulationModifier(refinement_modifier);


        simulator.Solve();

    }

};


#endif /*TESTPAPERSPEEDSIMULATION_HPP_*/
