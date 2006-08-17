/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/







# include <iostream>
# include <cstdlib>
# include <string>
# include <cstdio>


# include "NonSmoothLawFactory.h"


using namespace std;




namespace NonSmoothLawFactory
{


Registry& Registry :: get()
{



  static Registry instance ;
  return instance ;


}







void Registry :: add(const string& nom , object_creator creator)
{



  factory_map [ nom ] = creator ;

}




NonSmoothLaw* Registry :: instantiate(const std :: string& nom)
{

  std :: map < const std :: string , object_creator > :: iterator it = factory_map.find(nom) ;


  if (it != factory_map.end())         // found a factory?
  {

    return (it-> second)() ;      // run our factory

  }
  return 0 ;    // no factory found

}









Registration :: Registration(const string& nom ,  object_creator creator)
{


  cout << "Abonnement de " << nom << endl << endl ;

  Registry  ::get().add(nom,  creator) ;


}



}  // fin de namespace NonSmoothLawfactory








