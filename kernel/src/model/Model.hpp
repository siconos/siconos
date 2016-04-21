/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

/*! \file Model.hpp
  \brief The general Siconos Model
*/

#ifndef MODEL_H
#define MODEL_H

#include "SiconosConst.hpp"

#include "SiconosPointers.hpp"
#include "SiconosFwd.hpp"
#include "SiconosSerialization.hpp"

#include <string>


#ifndef __GCCXML__

#ifdef __GNUC__
#define DEPRECATED_MODEL_API(func) func __attribute__ ((deprecated ("This constructor of function is deprecrated and will be removed in the next major Siconos release! Use initialize() and setSimulation() !")))
#elif defined(_MSC_VER)
#define DEPRECATED_MODEL_API(func) __declspec(deprecated("This constructor will be removed in the next major Siconos release and does not work with MSVC 2013 !  Use initialize() and setSimulation() !")) func
#else
#define DEPRECATED_MODEL_API(func) func
#endif

#else
#define DEPRECATED_MODEL_API(func) func
#endif


/** \class Model
 * \brief  Model: object that links the NonSmoothDynamicalSystem with a
 * Simulation.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
 */
class Model : public std11::enable_shared_from_this<Model>
{
private:
  /** current time of the simulation 
      Warning FP : it corresponds to the time 
      at the end of the integration step. 
      It means that _t corresponds to tkp1 of the
      simulation or nextTime().
   */
  double _t;

  /** initial time of the simulation */
  double _t0;

  /** final time of the simulation */
  double _T;

  /** The simulation to solve the NonSmoothDynamicalSystem */
  SP::Simulation _strat;

  /** The NonSmoothDynamicalSystem of the simulation */
  SP::NonSmoothDynamicalSystem _nsds;

  /** information concerning the Model */
  std::string _title, _author, _description, _date;

  /** assignment operator => forbidden
      \return Model&
   */
  Model& operator=(const Model&);


  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(Model);

public:

  /** default constructor
   */
  Model();

     /** create the Model from a set of data
   *  \param t0 : the value for t0
   *  \param T : the value for T (optional parameter)
   *  \param title the title of the Model (optional parameter)
   *  \param author : the author of the Model (optional parameter)
   *  \param description : the description of the Model (optional
   *                  parameter)
   *  \param date : the date of the Model (optional parameter)
   */
  Model(double t0, double T = -1, const std::string& title = "none",
        const std::string& author = "nobody", const std::string& description = "none",
        const std::string& date = "none");

  /** destructor
   */
  ~Model();

  // --- GETTERS/SETTERS

  /** get the current time
   *  \return a double
   */
  inline double currentTime() const
  {
    return _t;
  }

  /** set the current time
   *  \param newValue the new time
   */
  inline void setCurrentTime(double newValue)
  {
    _t = newValue;
  }

  /** get initial time
   *  \return a double
   */
  inline double t0() const
  {
    return _t0;
  }

  /** set initial time of the time discretisation
   *  \param newT0
   */
  inline void sett0(double newT0)
  {
    _t0 = newT0;
  };

  /** get final time
   *  \return a double
   */
  inline double finalT() const
  {
    return _T;
  }

  /** set final time
   *  \param newValue the new final time for the Simulatiom
   */
  void setT(double newValue);

  /** get the Simulation of the Model
   *  \return a pointer on Simulation
   */
  inline SP::Simulation simulation() const
  {
    return _strat;
  }

  /** set the Simulation of the Model
   *  \return a pointer on Simulation
   */
  void setSimulation(SP::Simulation);

  /** get the NonSmoothDynamicalSystem of the Model
   *  \return a pointer on NonSmoothDynamicalSystem
   */
  inline SP::NonSmoothDynamicalSystem nonSmoothDynamicalSystem() const
  {
    return _nsds;
  }

  /** set the NonSmoothDynamicalSystem of the Model
   *  \param newPtr a pointer on NonSmoothDynamicalSystem
   */
  void setNonSmoothDynamicalSystemPtr(SP::NonSmoothDynamicalSystem newPtr);

  /** get the title of the simulation
   *  \return std::string : the title
   */
  inline const std::string  title() const
  {
    return _title;
  }

  /** set the title of the simulation
   *  \param s : the title
   */
  inline void setTitle(const std::string & s)
  {
    _title = s;
  }

  /** get the author of the simulation
   *  \return std::string : the author
   */
  inline const std::string  author() const
  {
    return _author;
  }

  /** set the author of the simulation
   *  \param s std::string : the author
   */
  inline void setAuthor(const std::string & s)
  {
    _author = s;
  }

  /** allows to get the description of the simulation
   *  \return std::string : the description
   */
  inline const std::string  description() const
  {
    return _description;
  }

  /** set the author of the simulation
   *  \param s std::string : the author
   */
  inline void setDescription(const std::string & s)
  {
    _description = s;
  }

  /** allows to get the date of the simulation
   *  \return std::string : the date
   */
  inline const std::string  date() const
  {
    return _date;
  }

  /** set the date of the simulation
   *  \param s std::string : the date
   */
  inline void setDate(const std::string & s)
  {
    _date = s;
  }

  /** Complete initialization of the model (NonSmoothDynamicalSystem,
      Simulation)
   */
  void initialize();

  DEPRECATED_MODEL_API(void initialize(SP::Simulation sim);)
  
   /** display the data of the Model
      \return void
   */
  void display() const ;
};

#endif // MODEL_H

