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









#ifndef _NONSMOOTHLAWFACTORY_H_        // This is supposed to avoid repeated inclusions
#define _NONSMOOTHLAWFACTORY_H_        // of the header file




# include <string>
# include <map>



# include "NonSmoothLaw.h"





class NonSmoothLaw ;





namespace NonSmoothLawFactory
{





template <  class SubType > NonSmoothLaw* factory()
{

  return new SubType ;

}


typedef NonSmoothLaw* (*object_creator)(void) ;      // base_creator est un pointeur sur une fonction renvoyant un pointeur de type Character*


class Registry
{


private :


  std :: map < const std :: string , object_creator > factory_map ;


public :

  static Registry& get() ;

  void add(const std :: string& nom , object_creator) ;

  NonSmoothLaw* instantiate(const std :: string& nom) ;


} ;


/******************************************************************

Classe permettant l'abonnement

*******************************************************************/

class Registration
{


public :

  Registration(const std :: string& nom , object_creator) ;


} ;






# define AUTO_REGISTER_NONSMOOTHLAW( class_name , class_type ) Registration  _registration_ ## class_type( class_name , &factory<class_type>);


// ATTENTION: ne pas mettre d'espace entre "AUTO_REGISTER_CHARACTER" et la parenthese ouvrante qui suit, sinon ca ne compile pas.



}  // fin de namespace CharacterFactory

#endif // ifndef _NONSMOOTHLAWFACTORY_H_













