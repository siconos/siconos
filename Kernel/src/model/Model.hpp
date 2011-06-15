/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

/*! \file Model.h
  \brief The general Siconos Model
*/

#ifndef MODEL_H
#define MODEL_H

#include "SiconosConst.hpp"
#include "Tools.hpp"
#include "SiconosPointers.hpp"
#include "InteractionsSet.hpp"

#include "SiconosSerialization.hpp"

class Simulation;
DEFINE_SPTR(SiconosModelXML);

/** Model: object that links the NonSmoothDynamicalSystem with a
 * Simulation.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
 */
class Model : public boost::enable_shared_from_this<Model>
{
private:
  /** current time of the simulation */
  double _t;

  /** initial time of the simulation */
  double _t0;

  /** final time of the simulation */
  double _T;

  /** The simulation to solve the NonSmoothDynamicalSystem */
  SP::Simulation _strat;

  /** The NonSmoothDynamicalSystem of the simulation */
  SP::NonSmoothDynamicalSystem _nsds;

  /** XML object linked to the Model */
  SP::SiconosModelXML _modelxml;

  /** information concerning the Model */
  std::string _title, _author, _description, _date, _xmlSchema;

  /** default constructor
   */
  Model();

  /** Copy constructor => private, no copy nor pass-by value for Model
   */
  Model(const Model&) {};

  /** assignment operator => forbidden
   */
  Model& operator=(const Model&);


  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(Model);

public:

  /** create the Model from an xml file
   *  \param char * : the input XML file (optional parameter)
   *  \exception RuntimeException
   */
  Model(const std::string& xmlFile);

  /** create the Model from a set of data
   *  \param double : the value for t0
   *  \param double : the value for T (optional parameter)
   *  \param string : the title of the Model (optional parameter)
   *  \param string : the author of the Model (optional parameter)
   *  \param string : the description of the Model (optional
   *                  parameter)
   *  \param string : the date of the Model (optional parameter)
   *  \param string : the xml schema of the Model (optional parameter)
   */
  Model(double, double = -1, const std::string& = "none",
        const std::string& = "nobody", const std::string& = "none",
        const std::string& = "none", const std::string& = "none");

  /** build the model from init/final times and a list of DS and
   *   Interactions
   *   \param double : the value for t0
   *   \param double : the value for T; if you do not want to set the
   *   final time, set T = -1.
   *   \param a list of DynamicalSystems
   *   \param a list of Interactions
   */
  Model(double, double, DynamicalSystemsSet&, InteractionsSet&);

  /** destructor
   */
  ~Model() {};

  // --- GETTERS/SETTERS

  /** get the current time
   *  \return a double
   */
  inline double currentTime() const
  {
    return _t;
  }

  /** set the current time
   *  \param a double
   */
  inline void setCurrentTime(const double& newValue)
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
   *  \param a double
   */
  inline void sett0(const double& newT0)
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
   *  \param a double
   */
  inline void setT(const double& newValue)
  {
    _T = newValue;
  }

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
  void setSimulationPtr(SP::Simulation);

  /** get the NonSmoothDynamicalSystem of the Model
   *  \return a pointer on NonSmoothDynamicalSystem
   */
  inline SP::NonSmoothDynamicalSystem nonSmoothDynamicalSystem() const
  {
    return _nsds;
  }

  /** set the NonSmoothDynamicalSystem of the Model
   *  \param a pointer on NonSmoothDynamicalSystem
   */
  void setNonSmoothDynamicalSystemPtr(SP::NonSmoothDynamicalSystem newPtr);

  /** get the SiconosModelXML of the Model
   *  \return a pointer on SiconosModelXML
   */
  inline SP::SiconosModelXML siconosModelXML() const
  {
    return _modelxml;
  }

  /** set the SiconosModelXML of the Model
   *  \param a pointer on SiconosModelXML
   */
  void setSiconosModelXMLPtr(SP::SiconosModelXML newPtr);

  /** get the title of the simulation
   *  \return string : the title
   */
  inline const std::string  title() const
  {
    return _title;
  }

  /** set the title of the simulation
   *  \param string : the title
   */
  inline void setTitle(const std::string & s)
  {
    _title = s;
  }

  /** get the author of the simulation
   *  \return string : the author
   */
  inline const std::string  author() const
  {
    return _author;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setAuthor(const std::string & s)
  {
    _author = s;
  }

  /** allows to get the description of the simulation
   *  \return string : the description
   */
  inline const std::string  description() const
  {
    return _description;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setDescription(const std::string & s)
  {
    _description = s;
  }

  /** allows to get the date of the simulation
   *  \return string : the date
   */
  inline const std::string  date() const
  {
    return _date;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setDate(const std::string & s)
  {
    _date = s;
  }

  /** allows to get the xmlSchema of the simulation
   *  \return string : the xmlSchema
   */
  inline const std::string  getXmlSchema() const
  {
    return _xmlSchema;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setXmlSchema(const std::string & s)
  {
    _xmlSchema = s;
  }

  /** Complete initialization of the model (NonSmoothDynamicalSystem,
      Simulation)
      \param a smart pointer to simulation (option, default = empty)
   */
  void initialize(SP::Simulation = SP::Simulation());

  // --- XML related functions ---

  /** saves into output file the data of the system
   *  \param char* : the data file which must be written
   */
  void saveToXMLFile(char*);

  /** saves into the DOM tree all the data of the system
   */
  void saveToDOMTree();

  /** copy the data of the plateform to the XML DOM tree
   *  \exception RuntimeException
   */
  void savePlatformToXML();

  /** check if the DOM tree respect the XML schema
   *  return bool : true if the DOM tree respect the XML schema
   *  \exception RuntimeException
   */
  bool checkXMLDOMTree();

  /** check if the XML objects of XML managment exist
   *  \exception RuntimeException
   */
  void checkXMLPlatform();

  /** check if the Model is complete. That's to say if the objects of
      the platform are coherent and if data of the XML are coherent
   *  \exception RuntimeException
   */
  void checkModelCoherency();

  /** checks if the xmlFile given respects the xmlSchema given
   *  \param string : the xml input file to check
   *  \param string : the xml schema
   *  \return int : 1 if the xml file respects the schema
   */
  int xmlSchemaValidated(std::string  xmlFile, std::string  xmlSchema = "");

  /** display the data of the Model
   */
  void display() const ;

};

#endif // MODEL_H

