/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

/*!\file ControlSimulation_impl.hpp
 * \brief functions related to the simulation involving control
 */

#ifndef ControlSimulation_impl_hpp
#define ControlSimulation_impl_hpp

#include <cstddef>
#include <string>
#include <utility>

#include "SiconosKernel.hpp"  // IWYU pragma: keep

#define TO_STR(x) std::to_string(x)

static inline std::pair<size_t, std::string> getNumberOfStates(
    siconos::graphs::DynamicalSystemsGraph& DSG0, siconos::graphs::InteractionsGraph& IG0) {
  std::string legend;
  siconos::graphs::DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  size_t nb = 0;
  int counter = 0;
  for (std::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi) {
    auto& x = *DSG0.bundle(*dsvi)->x();
    nb += x.size();

    std::string nameDS;
    if (DSG0.name.hasKey(*dsvi)) {
      nameDS = DSG0.name[*dsvi];
    } else {
      nameDS = "unknownDS" + TO_STR(counter);
      ++counter;
    }

    for (siconos::algebra::Index i = 0; i < x.size(); ++i) {
      legend.append(" " + nameDS + "_" + TO_STR(i));
    }

    if (DSG0.u.hasKey(*dsvi)) {
      siconos::algebra::Index sizeU = DSG0.u[*dsvi]->size();
      nb += sizeU;
      for (siconos::algebra::Index i = 0; i < sizeU; ++i) {
        legend.append(" " + nameDS + "_u_" + TO_STR(i));
      }
    }

    if (DSG0.eVector.hasKey(*dsvi)) {
      siconos::algebra::Index sizeE = DSG0.eVector[*dsvi]->size();
      for (siconos::algebra::Index i = 0; i < sizeE; ++i) {
        legend.append(" " + nameDS + "_e_" + TO_STR(i));
      }
      nb += DSG0.eVector[*dsvi]->size();
    }
  }

  siconos::graphs::InteractionsGraph::VIterator ivi, ivdend;
  counter = 0;
  for (std::tie(ivi, ivdend) = IG0.vertices(); ivi != ivdend; ++ivi) {
    std::string nameInter;
    if (IG0.name.hasKey(*ivi)) {
      nameInter = IG0.name[*ivi];
    } else {
      nameInter = "unknownInteraction" + TO_STR(counter);
      ++counter;
    }
    auto& y = *IG0.bundle(*ivi)->y(0);
    nb += y.size();
    for (siconos::algebra::Index i = 0; i < y.size(); ++i) {
      legend.append(" " + nameInter + "_y_" + TO_STR(i));
    }

    auto& lambda = *IG0.bundle(*ivi)->lambda(0);
    nb += lambda.size();
    for (siconos::algebra::Index i = 0; i < lambda.size(); ++i) {
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
 * \return the index of the last written column
 */
static inline size_t storeAllStates(size_t indx, size_t startColumn,
                                    siconos::graphs::DynamicalSystemsGraph& DSG0,
                                    siconos::graphs::InteractionsGraph& IG0,
                                    siconos::algebra::SiconosMatrix& data) {
  siconos::graphs::DynamicalSystemsGraph::VIterator dsvi, dsvdend;
  size_t column = startColumn;
  for (std::tie(dsvi, dsvdend) = DSG0.vertices(); dsvi != dsvdend; ++dsvi) {
    size_t i = column;
    auto& x = *DSG0.bundle(*dsvi)->x();
    for (siconos::algebra::Index j = 0; j < x.size(); ++i, ++j) {
      data(indx, i) = x(j);
    }
    column += x.size();

    if (DSG0.u.hasKey(*dsvi)) {
      auto& u = *DSG0.u[*dsvi];
      for (siconos::algebra::Index j = 0; j < u.size(); ++i, ++j) {
        data(indx, i) = u(j);
      }
      column += u.size();
    }

    if (DSG0.eVector.hasKey(*dsvi)) {
      auto& e = *DSG0.eVector[*dsvi];
      for (siconos::algebra::Index j = 0; j < e.size(); ++i, ++j) {
        data(indx, i) = e(j);
      }
      column += e.size();
    }
  }

  siconos::graphs::InteractionsGraph::VIterator ivi, ivdend;
  for (std::tie(ivi, ivdend) = IG0.vertices(); ivi != ivdend; ++ivi) {
    size_t i = column;
    auto& y = *IG0.bundle(*ivi)->y(0);
    for (siconos::algebra::Index j = 0; j < y.size(); ++i, ++j) {
      data(indx, i) = y(j);
    }
    column += y.size();

    auto& lambda = *IG0.bundle(*ivi)->lambda(0);
    for (siconos::algebra::Index j = 0; j < lambda.size(); ++i, ++j) {
      data(indx, i) = lambda(j);
    }
    column += lambda.size();
  }

  return column;
}

#endif
