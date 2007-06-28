/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

/*! \file*/

#ifndef MODEL_H
#define MODEL_H

#include "SiconosConst.h"
#include "Tools.h"

class NonSmoothDynamicalSystem;
class Simulation;
class SiconosModelXML;
class TimeDiscretisation;

/** Objects that links the NonSmoothDynamicalSystem with its Simulation.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.0.
 *  \date (Creation) Apr 26, 2004
 *
 *
 */
class Model
{
private:
  /** current time of the simulation */
  double t;

  /** initial time of the simulation */
  double t0;

  /** final time of the simulation */
  double T;

  /** The simulation to solve the NonSmoothDynamicalSystem */
  Simulation *strat;

  /** The NonSmoothDynamicalSystem of the simulation */
  NonSmoothDynamicalSystem * nsds;

  /** XML object linked to the Model */
  SiconosModelXML *modelxml;

  /** information concerning the Model */
  std::string  title, author, description, date, xmlSchema;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  BoolMap isAllocatedIn;

  /** default constructor
   */
  Model();

  /** Copy constructor => private, no copy nor pass-by value for Model
   */
  Model(const Model&);

public:

  /** create the Model from an xml file
   *  \param char * : the input XML file (optional parameter)
   *  \exception RuntimeException
   */
  Model(char *xmlFile);

  /** create the Model from a set of data
   *  \param double : the value for t0
   *  \param double : the value for T (optional parameter)
   *  \param string : the title of the Model (optional parameter)
   *  \param string : the author of the Model (optional parameter)
   *  \param string : the description of the Model (optional parameter)
   *  \param string : the date of the Model (optional parameter)
   *  \param string : the xml schema of the Model (optional parameter)
   *  \exception RuntimeException
   */
  Model(double, double = -1, const std::string& = "none", const std::string& = "nobody",
        const std::string& = "none", const std::string& = "none", const std::string& = "none");

  /** destructor
   */
  ~Model();

  // --- GETTERS/SETTERS

  /** get the current time
   *  \return a double
   */
  inline const double getCurrentT() const
  {
    return t;
  }

  /** set the current time
   *  \param a double
   */
  inline void setCurrentT(const double& newValue)
  {
    t = newValue;
  }

  /** get initial time
   *  \return a double
   */
  inline const double getT0() const
  {
    return t0;
  }

  /** set initial time of the time discretisation
   *  \param a double
   */
  inline void setT0(const double& newT0)
  {
    t0 = newT0;
  };

  /** get final time
   *  \return a double
   */
  inline const double getFinalT() const
  {
    return T;
  }

  /** set final time
   *  \param a double
   */
  inline void setFinalT(const double& newValue)
  {
    T = newValue;
  }

  /** get the Simulation of the Model
   *  \return a pointer on Simulation
   */
  inline Simulation* getSimulationPtr() const
  {
    return strat;
  }

  /** set the Simulation of the Model
   *  \return a pointer on Simulation
   */
  void setSimulationPtr(Simulation*);

  /** get the NonSmoothDynamicalSystem of the Model
   *  \return a pointer on NonSmoothDynamicalSystem
   */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
  {
    return nsds;
  }

  /** set the NonSmoothDynamicalSystem of the Model
   *  \param a pointer on NonSmoothDynamicalSystem
   */
  void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newPtr);

  /** get the SiconosModelXML of the Model
   *  \return a pointer on SiconosModelXML
   */
  inline SiconosModelXML* getSiconosModelXMLPtr() const
  {
    return modelxml;
  }

  /** set the SiconosModelXML of the Model
   *  \param a pointer on SiconosModelXML
   */
  void setSiconosModelXMLPtr(SiconosModelXML *newPtr);

  /** get the title of the simulation
   *  \return string : the title
   */
  inline const std::string  getTitle() const
  {
    return title;
  }

  /** set the title of the simulation
   *  \param string : the title
   */
  inline void setTitle(const std::string & s)
  {
    title = s;
  }

  /** get the author of the simulation
   *  \return string : the author
   */
  inline const std::string  getAuthor() const
  {
    return author;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setAuthor(const std::string & s)
  {
    author = s;
  }

  /** allows to get the description of the simulation
   *  \return string : the description
   */
  inline const std::string  getDescription() const
  {
    return description;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setDescription(const std::string & s)
  {
    description = s;
  }

  /** allows to get the date of the simulation
   *  \return string : the date
   */
  inline const std::string  getDate() const
  {
    return date;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setDate(const std::string & s)
  {
    date = s;
  }

  /** allows to get the xmlSchema of the simulation
   *  \return string : the xmlSchema
   */
  inline const std::string  getXmlSchema() const
  {
    return xmlSchema;
  }

  /** set the author of the simulation
   *  \param string : the author
   */
  inline void setXmlSchema(const std::string & s)
  {
    xmlSchema = s;
  }

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

  /** check if the Model is complete. That's to say if the objects of the platform are coherent and if data of the XML are coherent
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

