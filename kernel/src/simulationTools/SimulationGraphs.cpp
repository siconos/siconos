/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2017 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "SimulationGraphs.hpp"
#include "DynamicalSystem.hpp"
#include "SimpleMatrix.hpp"
#include "Interaction.hpp"

unsigned int InteractionsGraph::update_positions()
{
  unsigned dim = 0;
  VIterator vd, vdend;
  for (std11::tie(vd, vdend) = vertices(); vd != vdend; ++vd)
    {
      properties(*vd).absolute_position = dim;
      dim += (bundle(*vd)->dimension());
    }
  return dim;
}


void InteractionsGraph::display_edges_blocks()
{
  InteractionsGraph::EIterator ei, eiend;
  unsigned int isource, itarget;
  unsigned int source_size, target_size;
  bool compute_upper;
  InteractionsGraph::VDescriptor vd_source, vd_target;
  SP::Interaction intsource, intarget;

  // Loop over all edges (dynamical systems!) in active interactions set
  for (std11::tie(ei, eiend) = edges(); ei != eiend; ++ei)
    {
      vd_source = source(*ei);
      vd_target = target(*ei);
      intsource = bundle(vd_source);
      intarget = bundle(vd_target);
      // ed_parent = parentSet.descriptor(bundle(*ei));
      // their sizes, 
      // and their indices.
      isource = index(vd_source);
      itarget = index(vd_target);
      compute_upper = (itarget > isource);
      //std::cout << "DEAL WITH EDGE " << index(*ei) << ", ds number : " << bundle(*ei)->number()<< std::endl;
      //std::cout << "BEtween (index/inter number)" << isource << "/" << source->number() << " and " << itarget << "/" << target->number() << std::endl;
      if (compute_upper) // upper block
	{
	  if ( properties(*ei).upper_block)
	    {
	      std::cout << "UPPER, DEAL WITH EDGE " << index(*ei) << ", ds number : " << bundle(*ei)->number()<< std::endl;
	      std::cout << "BEtween (index/inter number)" << isource << "/" << intsource->number() << " and " << itarget << "/" << intarget->number() << std::endl;
	      //std::cout << "ED lin upper block allocation " << isource << " " << itarget << " " << bundle(*ei)->number() << std::endl;
	      properties(*ei).upper_block->display();
	      properties(*ei).lower_block->display();
	      std::cout << " Pointer upper value " << properties(*ei).upper_block << std::endl;
	    }
	}
      else // lower block
	{
	  if ( properties(*ei).lower_block)
	    {
	      std::cout << "LOWER, DEAL WITH EDGE " << index(*ei) << ", ds number : " << bundle(*ei)->number()<< std::endl;
	      std::cout << "BEtween (index/inter number)" << isource << "/" << intsource->number() << " and " << itarget << "/" << intarget->number() << std::endl;
	      properties(*ei).lower_block->display();
	      properties(*ei).upper_block->display();
	      std::cout << " Pointer lower value " << properties(*ei).lower_block << std::endl;
	    }
	}
    }
    
}

