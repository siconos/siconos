#ifndef MODEL_H
#define MODEL_H

#include "SiconosModelXML.h"
#include "NonSmoothDynamicalSystem.h"
#include "Strategy.h"

#include "check.h"
#include "SiconosConst.h"

#include <iostream>
#include <vector>

extern unsigned int MATRIX_MAX_SIZE;
extern unsigned int VECTOR_MAX_SIZE;
extern std::string  FILE_STORAGE;
extern std::string  XML_SCHEMA;

class NonSmoothDynamicalSystem;
class Strategy;
class SiconosModelXML;

/** \class Model
 *  \brief The Model regroups the high level functionnalities of the platform
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 26, 2004
 *
 *
 */
class Model
{
public:

  /** \fn Model(char *xmlFile)
   *  \brief create the Model from an xml file
   *  \param char * : the input XML file (optional parameter)
   *  \exception RuntimeException
   */
  Model(char *xmlFile);

  /** \fn Model(double t0, double T, string title, string author, string description, string date, string xmlSchema)
   *  \brief create the Model from a set of data
   *  \param double : the value for t0
   *  \param double : the value for T (optional parameter)
   *  \param string : the title of the Model (optional parameter)
   *  \param string : the author of the Model (optional parameter)
   *  \param string : the description of the Model (optional parameter)
   *  \param string : the date of the Model (optional parameter)
   *  \param string : the xml schema of the Model (optional parameter)
   *  \exception RuntimeException
   */
  Model(double t0, double T = -1, std::string  title = "none", std::string  author = "nobody",
        std::string  description = "none", std::string  date = "none",
        std::string  xmlSchema = "none");

  /** \fn Model(double t0, double T, string title, string author, string description, string date, string xmlSchema)
   *  \brief create the Model from a set of data
   *  \param double : the value for t0
   *  \param string : the type of strategy to be used
   *  \param double : the value for T (optional parameter)
   *  \param string : the title of the Model (optional parameter)
   *  \param string : the author of the Model (optional parameter)
   *  \param string : the description of the Model (optional parameter)
   *  \param string : the date of the Model (optional parameter)
   *  \param string : the xml schema of the Model (optional parameter)
   *  \exception RuntimeException
   */
  /*Model(double t0, string strategy_name, double T = -1, string title="none", string author="nobody",
    string description="none", string date="none",
    string xmlSchema="none");
  */
  /** \fn Model(double t0, double T, string title, string author, string description, string date, string xmlSchema)
   *  \brief create the Model from a set of data
   *  \param double : the value for t0
   *  \param strategy : the predefined strategy to use
   *  \param double : the value for T (optional parameter)
   *  \param string : the title of the Model (optional parameter)
   *  \param string : the author of the Model (optional parameter)
   *  \param string : the description of the Model (optional parameter)
   *  \param string : the date of the Model (optional parameter)
   *  \param string : the xml schema of the Model (optional parameter)
   *  \exception RuntimeException
   */
  /*Model(double t0, string strategy_name, double T = -1, string title="none", string author="nobody",
    string description="none", string date="none",
    string xmlSchema="none");
  */
  ~Model();

  // --- GETTERS/SETTERS

  /** \fn inline const double getCurrentT()
   *  \brief get the current time
   *  \return a double
   */
  inline const double getCurrentT() const
  {
    return t;
  }

  /** \fn inline void setCurrentT(const double&)
   *  \brief set the current time
   *  \param a double
   */
  inline void setCurrentT(const double& newValue)
  {
    t = newValue;
  }

  /** \fn inline const double getT0()
   *  \brief get initial time
   *  \return a double
   */
  inline const double getT0() const
  {
    return t0;
  }

  /** \fn inline void setT0(const double&)
   *  \brief set initial time
   *  \param a double
   */
  inline void setT0(const double& newValue)
  {
    t0 = newValue;
  }

  /** \fn inline const double getFinalT()
   *  \brief get final time
   *  \return a double
   */
  inline const double getFinalT() const
  {
    return T;
  }

  /** \fn inline void setFinalT(const double&)
   *  \brief set final time
   *  \param a double
   */
  inline void setFinalT(const double& newValue)
  {
    T = newValue;
  }

  /** \fn inline Strategy* getStrategyPtr() const
   *  \brief get the Strategy of the Model
   *  \return a pointer on Strategy
   */
  inline Strategy* getStrategyPtr() const
  {
    return strat;
  }

  /** \fn void setStrategyPtr(Strategy *str)
   *  \brief set the Strategy of the Model
   *  \param a pointer on Strategy
   */
  void setStrategyPtr(Strategy * newPtr);


  //  void linkModelToStrategy(model *, strategy *)

  /** \fn inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
   *  \brief get the NonSmoothDynamicalSystem of the Model
   *  \return a pointer on NonSmoothDynamicalSystem
   */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
  {
    return nsds;
  }

  /** \fn void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newPtr)
   *  \brief set the NonSmoothDynamicalSystem of the Model
   *  \param a pointer on NonSmoothDynamicalSystem
   */
  void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newPtr);

  /** \fn inline SiconosModelXML* getSiconosModelXMLPtr() const
   *  \brief get the SiconosModelXML of the Model
   *  \return a pointer on SiconosModelXML
   */
  inline SiconosModelXML* getSiconosModelXMLPtr() const
  {
    return modelxml;
  }

  /** \fn void setSiconosModelXMLPtr(SiconosModelXML *newPtr)
   *  \brief set the SiconosModelXML of the Model
   *  \param a pointer on SiconosModelXML
   */
  void setSiconosModelXMLPtr(SiconosModelXML *newPtr);

  /** \fn inline string getTitle()
   *  \brief get the title of the simulation
   *  \return string : the title
   */
  inline const std::string  getTitle() const
  {
    return title;
  }

  /** \fn inline void setTitle(const string& s)
   *  \brief set the title of the simulation
   *  \param string : the title
   */
  inline void setTitle(const std::string & s)
  {
    title = s;
  }

  /** \fn inline string getAuthor()
   *  \brief get the author of the simulation
   *  \return string : the author
   */
  inline const std::string  getAuthor() const
  {
    return author;
  }

  /** \fn inline void setAuthor(const string& s)
   *  \brief set the author of the simulation
   *  \param string : the author
   */
  inline void setAuthor(const std::string & s)
  {
    author = s;
  }

  /** \fn inline const string getDescription()
   *  \brief allows to get the description of the simulation
   *  \return string : the description
   */
  inline const std::string  getDescription() const
  {
    return description;
  }

  /** \fn inline void setDescription(const string& s)
   *  \brief set the author of the simulation
   *  \param string : the author
   */
  inline void setDescription(const std::string & s)
  {
    description = s;
  }

  /** \fn inline const string getDate()
   *  \brief allows to get the date of the simulation
   *  \return string : the date
   */
  inline const std::string  getDate() const
  {
    return date;
  }

  /** \fn inline void setDate(const string& s)
   *  \brief set the author of the simulation
   *  \param string : the author
   */
  inline void setDate(const std::string & s)
  {
    date = s;
  }

  /** \fn inline const string getXmlSchema()
   *  \brief allows to get the xmlSchema of the simulation
   *  \return string : the xmlSchema
   */
  inline const std::string  getXmlSchema() const
  {
    return xmlSchema;
  }

  /** \fn inline void setXmlSchema(const string& s)
   *  \brief set the author of the simulation
   *  \param string : the author
   */
  inline void setXmlSchema(const std::string & s)
  {
    xmlSchema = s;
  }

  // --- XML related functions ---

  /** \fn saveToXMLFile(char*)
   *  \brief saves into output file the data of the system
   *  \param char* : the data file which must be written
   */
  void saveToXMLFile(char*);

  /** \fn saveToDOMTree()
   *  \brief saves into the DOM tree all the data of the system
   */
  void saveToDOMTree();

  /** \fn void savePlatformToXML()
   *  \brief copy the data of the plateform to the XML DOM tree
   *  \exception RuntimeException
   */
  void savePlatformToXML();

  /** \fn bool checkXMLDOMTree()
   *  \brief check if the DOM tree respect the XML schema
   *  return bool : true if the DOM tree respect the XML schema
   *  \exception RuntimeException
   */
  bool checkXMLDOMTree();

  /** \fn void checkXMLPlatform()
   *  \brief check if the XML objects of XML managment exist
   *  \exception RuntimeException
   */
  void checkXMLPlatform();

  /** \fn void checkModelCoherency()
   *  \brief check if the Model is complete. That's to say if the objects of the platform are coherent and if data of the XML are coherent
   *  \exception RuntimeException
   */
  void checkModelCoherency();

  /** \fn int xmlSchemaValidated(string xmlFile, string xmlSchema)
   *  \brief checks if the xmlFile given respects the xmlSchema given
   *  \param string : the xml input file to check
   *  \param string : the xml schema
   *  \return int : 1 if the xml file respects the schema
   */
  int xmlSchemaValidated(std::string  xmlFile, std::string  xmlSchema = "");

  // --- OTHER FUNCTIONS ---

  /** \fn bool isModelComplete(void)
   *  \brief determines if there are enough data to formalise and solve the problem
   *  \return tru if the Model is complete
   */
  //bool isModelComplete();

  /** \fn void runSimulation(void)
   *  \brief launches the simulation
   */
  //void runSimulation(void){};

  /** \fn void doOneStep(void)
   *  \brief makes the simulation go on in time
   */
  //void doOneStep(void);

  /** \fn void display()
   *  \brief display the data of the Model
   */
  void display() const ;

  /*******************************************************
   *
   * function to create the platform from a C++ programm
   *
   * \todo : interface methods of the Model
   *
   *//////////////////////////////////////////////////////

  /** \fn NonSmoothDynamicalSystem* createNonSmoothDynamicalSystem(bool bvp)
   *  \brief allows to create the NonSmoothDynamicalSystem of the Model
   *  \param bool : to define if the NonSmoothDynamicalSystem is BVP or not
   *  \return NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem created
   */
  NonSmoothDynamicalSystem* createNonSmoothDynamicalSystem(bool bvp = false);

  /** \fn Strategy* createStrategy(string type)
   *  \brief allows to create a Strategy (EventDriven or TimeStepping)
   *  \return Strategy* : the Strategy created
   */
  Strategy* createStrategy(std::string  type);

  /** \fn Strategy* createTimeStepping()
   *  \brief allows to create a Strategy : TimeStepping
   *  \return Strategy* : the Strategy created
   */
  Strategy* createTimeStepping();

  /** \fn Strategy* createTimeEventDriven()
   *  \brief allows to create a Strategy : EventDriven
   *  \return Strategy* : the Strategy created
   */
  Strategy* createTimeEventDriven();

  /** \fn void setMatrixMaxSize( int max )
   *  \brief allows to change the value MatrixMaxSize which determines the bound over between
   * external file storing and xml input/outpout file storing
   *  \param int :  the value to assign to MatrixMaxSize
   */
  void setMatrixMaxSize(const int& max);

  /** \fn void setVectorMaxSize( int max )
   *  \brief allows to change the value VectorMaxSize which determines the bound over between
   * external file storing and xml input/outpout file storing
   *  \param int :  the value to assign to VectorMaxSize
   */
  void setVectorMaxSize(const int& max);

  /** \fn void setFileStorage( string fs )
   *  \brief allows to change the value fileStorage which determines the format
   * of external file save (binary or ascii)
   *  \param string :  the value to assign to fileStorage
   */
  void setFileStorage(const std::string & fs);

private:
  /** \fn Model()
   *  \brief default constructor
   */
  Model();

  /** current time of the simulation */
  double t;

  /** initial time of the simulation */
  double t0;

  /** final time of the simulation */
  double T;

  /** The strategy to solve the NonSmoothDynamicalSystem */
  Strategy *strat;

  /** The NonSmoothDynamicalSystem of the simulation */
  NonSmoothDynamicalSystem * nsds;

  /** XML object linked to the Model */
  SiconosModelXML *modelxml;

  /** information concerning the Model */
  std::string  title, author, description, date, xmlSchema;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isNsdsAllocatedIn;
  bool isStrategyAllocatedIn;
  bool isModelXmlAllocatedIn;

};

#endif // MODEL_H

