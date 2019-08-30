/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "frictionContact_test_utils.h"

char ** data_collection()
{

  int n_data_1=150;

  char ** data_collection_1 = (char **)malloc(n_data_1*sizeof(char *));
  int n_data=0;


  data_collection_1[n_data++] = "./data/FC3D_Example1_SBM.dat";
  data_collection_1[n_data++] = "./data/Capsules-i122-1617.dat";
  data_collection_1[n_data++] = "./data/Confeti-ex13-4contact-Fc3D-SBM.dat";
  data_collection_1[n_data++] = "./data/GFC3D_TwoRods1-condensed.dat";
  data_collection_1[n_data++] = "./data/FC3D_Example1.dat";
  data_collection_1[n_data++] = "./data/Confeti-ex03-Fc3D-SBM.dat";         /*5*/
  data_collection_1[n_data++] = "./data/BoxesStack1-i100000-32.hdf5.dat";   
  data_collection_1[n_data++] = "./data/RockPile_tob1.dat";
  data_collection_1[n_data++] = "./data/KaplasTower-i1061-4.hdf5.dat";
  data_collection_1[n_data++] = "./data/OneObject-i100000-499.hdf5.dat";
  data_collection_1[n_data++] = "./data/Rover4396.dat";                    /*10*/
  data_collection_1[n_data++] = "---";




  return data_collection_1;
}
