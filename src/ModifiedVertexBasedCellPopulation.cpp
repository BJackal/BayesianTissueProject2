/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "ModifiedVertexBasedCellPopulation.hpp"
#include <boost/foreach.hpp>
#include "VertexMeshWriter.hpp"
#include "Warnings.hpp"
#include "ChasteSyscalls.hpp"
#include "IsNan.hpp"
#include "ShortAxisVertexBasedDivisionRule.hpp"

// Cell writers
#include "CellAgesWriter.hpp"
#include "CellAncestorWriter.hpp"
#include "CellProliferativePhasesWriter.hpp"
#include "CellProliferativeTypesWriter.hpp"
#include "CellVolumesWriter.hpp"

// Cell population writers
#include "CellMutationStatesCountWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "VertexT1SwapLocationsWriter.hpp"
#include "VertexT2SwapLocationsWriter.hpp"
#include "VertexT3SwapLocationsWriter.hpp"
#include "StepSizeException.hpp"

#include <iostream>
#include <fstream>

template<unsigned DIM>
ModifiedVertexBasedCellPopulation<DIM>::ModifiedVertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh,
                                          std::vector<CellPtr>& rCells,
                                          bool deleteMesh,
                                          bool validate,
                                          const std::vector<unsigned> locationIndices)
    : VertexBasedCellPopulation<DIM>(rMesh, rCells, deleteMesh, validate, locationIndices),
	  mRestrictVertexMovement(true)
{
}

template<unsigned DIM>
ModifiedVertexBasedCellPopulation<DIM>::ModifiedVertexBasedCellPopulation(MutableVertexMesh<DIM, DIM>& rMesh)
    : VertexBasedCellPopulation<DIM>(rMesh)
{
}

template<unsigned DIM>
ModifiedVertexBasedCellPopulation<DIM>::~ModifiedVertexBasedCellPopulation()
{
}

template<unsigned DIM>
void ModifiedVertexBasedCellPopulation<DIM>::SetRestrictVertexMovementBoolean(bool restrictVertexMovement)
{
    mRestrictVertexMovement = restrictVertexMovement;
}

template<unsigned DIM>
bool ModifiedVertexBasedCellPopulation<DIM>::GetRestrictVertexMovementBoolean()
{
	return mRestrictVertexMovement;
}

template<unsigned DIM>
void ModifiedVertexBasedCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<ModifiedVertexBasedCellPopulation>\n";
    *rParamsFile << "\t\t\t<RestrictVertexMovement>" << this->mRestrictVertexMovement << "</RestrictVertexMovement>\n";
    *rParamsFile << "\t\t</ModifiedVertexBasedCellPopulation>\n";

    // Call method on direct parent class
    VertexBasedCellPopulation<DIM>::OutputCellPopulationParameters(rParamsFile);
}

template<unsigned DIM>
void ModifiedVertexBasedCellPopulation<DIM>::CheckForStepSizeException(unsigned nodeIndex, c_vector<double,DIM>& rDisplacement, double dt)
{
    double length = norm_2(rDisplacement);

    if (this->mRestrictVertexMovement)
    {
        if (length > 0.5*this->mpMutableVertexMesh->GetCellRearrangementThreshold())
        {
            rDisplacement *= 0.5*this->mpMutableVertexMesh->GetCellRearrangementThreshold()/length;

            std::ostringstream message;
            message << "Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted ";
            message << "so the motion has been restricted. Use a smaller timestep to avoid these warnings.";

            double suggested_step = 0.95*dt*((0.5*this->mpMutableVertexMesh->GetCellRearrangementThreshold())/length);

            throw StepSizeException(suggested_step, message.str(), false);
        }
    }
}

template<unsigned DIM>
void ModifiedVertexBasedCellPopulation<DIM>::UpdateNodeLocations(double dt)
{
    // Iterate over all nodes associated with real cells to update their positions
    for (unsigned node_index=0; node_index<this->GetNumNodes(); node_index++)
    {
        // Get the damping constant for this node
        double damping_const = this->GetDampingConstant(node_index);

        // Compute the displacement of this node
        c_vector<double, DIM> displacement = dt*this->GetNode(node_index)->rGetAppliedForce()/damping_const;

        if(mRestrictVertexMovement)
        {
            /*
             * If the displacement of this node is greater than half the cell rearrangement threshold,
             * this could result in nodes moving into the interior of other elements, which should not
             * be possible. Therefore in this case we restrict the displacement of the node to the cell
             * rearrangement threshold and warn the user that a smaller timestep should be used. This
             * restriction ensures that vertex elements remain well defined (see #1376).
             */
            if (norm_2(displacement) > 0.5*this->rGetMesh().GetCellRearrangementThreshold())
            {
                WARN_ONCE_ONLY("Vertices are moving more than half the CellRearrangementThreshold. This could cause elements to become inverted so the motion has been restricted. Use a smaller timestep to avoid these warnings.");
                displacement *= 0.5*this->rGetMesh().GetCellRearrangementThreshold()/norm_2(displacement);
            }
        }

        // Get new node location
        c_vector<double, DIM> new_node_location = this->GetNode(node_index)->rGetLocation() + displacement;

        for (unsigned i=0; i<DIM; i++)
        {
            assert(!std::isnan(new_node_location(i)));
        }

        // Create ChastePoint for new node location
        ChastePoint<DIM> new_point(new_node_location);

        // Move the node
        this->SetNode(node_index, new_point);
    }
}

// Explicit instantiation
template class ModifiedVertexBasedCellPopulation<1>;
template class ModifiedVertexBasedCellPopulation<2>;
template class ModifiedVertexBasedCellPopulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ModifiedVertexBasedCellPopulation)
