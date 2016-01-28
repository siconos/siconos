/* Siconos-Kernel, Copyright INRIA 2005-2014
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*!\file ControlSimulation_impl.hpp
 * \brief functions related to the simulation involving control
 */

#ifndef ControlSimulation_impl_hpp
#define ControlSimulation_impl_hpp

#include <utility>

#include "SimulationTypeDef.hpp"

#include <SiconosConfig.h>
#if defined(SICONOS_STD_TO_STRING) && !defined(SICONOS_USE_BOOST_FOR_CXX11)
#define TO_STR(x) std::to_string(x)
#else
#include <boost/lexical_cast.hpp>
#define TO_STR(x) boost::lexical_cast<std::string>(x)
#endif


static inline std::pair<unsigned, std::string> getNumberOfStates(DynamicalSystemsGraph& DSG0, InteractionsGraph& IG0)
{
  std::string legend;
  DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  unsigned nb = 0;
  unsigned counter = 0;
  for (std11::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi)
  {
    SiconosVector& x = *DSG0.bundle(*dsvi)->x();
    nb += x.size();
    std::string nameDS;
    if (DSG0.name.hasKey(*dsvi))
    {
      nameDS = DSG0.name[*dsvi];
    }
    else
    {
      nameDS = "unknownDS" + TO_STR(counter);
      ++counter;
    }

    for (unsigned i = 0; i < x.size(); ++i)
    {
      legend.append(" " + nameDS + "_" + TO_STR(i));
    }



    if (DSG0.u.hasKey(*dsvi))
    {
      unsigned sizeU = DSG0.u[*dsvi]->size();
      nb += sizeU;
      for (unsigned i = 0; i < sizeU; ++i)
      {
        legend.append(" " + nameDS + "_u_" + TO_STR(i));
      }
    }

    if (DSG0.e.hasKey(*dsvi))
    {
      unsigned sizeE = DSG0.e[*dsvi]->size();
      for (unsigned i = 0; i < sizeE; ++i)
      {
        legend.append(" " + nameDS + "_e_" + TO_STR(i));
      }
      nb += DSG0.e[*dsvi]->size();
    }
  }

  InteractionsGraph::VIterator ivi, ivdend;
  counter = 0;
  for (std11::tie(ivi, ivdend) = IG0.vertices(); ivi != ivdend; ++ivi)
  {
    std::string nameInter;
    if (IG0.name.hasKey(*ivi))
    {
      nameInter = IG0.name[*ivi];
    }
    else
    {
      nameInter = "unknownInteraction" + TO_STR(counter);
      ++counter;
    }
    SiconosVector& y = *IG0.bundle(*ivi)->y(0);
    nb += y.size();
    for (unsigned i = 0; i < y.size(); ++i)
    {
      legend.append(" " + nameInter + "_y_" + TO_STR(i));
    }

    SiconosVector& lambda = *IG0.bundle(*ivi)->lambda(0);
    nb += lambda.size();
    for (unsigned i = 0; i < lambda.size(); ++i)
    {
      legend.append(" " + nameInter + "_lambda_" + TO_STR(i));
    }
  }

  return std::make_pair(nb, legend);
}

/** store all the states of the graph in a matrix
 * \param indx row index in the matrix
 * \param startColumn the starting column
 * \param DSG0 the graph of DynamicalSystem
 * \param IG0 the graph of Interaction
 * \param data the matrix where to save the data
 * \return the last written column
 */
static inline unsigned storeAllStates(unsigned indx, unsigned startColumn, DynamicalSystemsGraph& DSG0, InteractionsGraph& IG0, SimpleMatrix& data)
{
  DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  unsigned column = startColumn;
  for (std11::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi)
  {
    unsigned i = column;
    SiconosVector& x = *DSG0.bundle(*dsvi)->x();
    for (unsigned j = 0; j < x.size(); ++i, ++j)
    {
      data(indx, i) = x(j);
    }
    column += x.size();

    if (DSG0.u.hasKey(*dsvi))
    {
      SiconosVector& u = *DSG0.u[*dsvi];
      for (unsigned j = 0; j < u.size(); ++i, ++j)
      {
        data(indx, i) = u(j);
      }
      column += u.size();
    }

    if (DSG0.e.hasKey(*dsvi))
    {
      SiconosVector& e = *DSG0.e[*dsvi];
      for (unsigned j = 0; j < e.size(); ++i, ++j)
      {
        data(indx, i) = e(j);
      }
      column += e.size();
    }

  }

  InteractionsGraph::VIterator ivi, ivdend;
  for (std11::tie(ivi, ivdend) = IG0.vertices(); ivi != ivdend; ++ivi)
  {
    unsigned i = column;
    SiconosVector& y = *IG0.bundle(*ivi)->y(0);
    for (unsigned j = 0; j < y.size(); ++i, ++j)
    {
      data(indx, i) = y(j);
    }
    column += y.size();

    SiconosVector& lambda = *IG0.bundle(*ivi)->lambda(0);
    for (unsigned j = 0; j < lambda.size(); ++i, ++j)
    {
      data(indx, i) = lambda(j);
    }
    column += lambda.size();
  }

  return column;
}

#endif
