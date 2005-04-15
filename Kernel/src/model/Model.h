#ifndef MODEL_H
#define MODEL_H


#include <iostream>
#include <vector>
#include "SiconosModelXML.h"
#include "NonSmoothDynamicalSystem.h"
#include "Strategy.h"

#include "SiconosConst.h"
//#include "KernelDefaultConfig.h"


//using namespace std;

extern int MATRIX_MAX_SIZE;
extern int VECTOR_MAX_SIZE;
extern string FILE_STORAGE;
extern string XML_SCHEMA;

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
  /** \fn Model()
   *  \brief constructor by default. Put pointers to NULL.
   */
  Model();

  //  /** \fn Model(char *xmlFile, float t=-1, float t0=-1, float T=-1, NonSmoothDynamicalSystem* nsds=NULL, Strategy* strategy=NULL)
  //   *  \brief constructor which creates the XML structure by reading the xml file given in parameter
  //   *  \param char * : the input XML file
  //   *  \param float : the value for t (optional parameter)
  //   *  \param float : the value for t0 (optional parameter)
  //   *  \param float : the value for T (optional parameter)
  //   *  \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem linked to the Model (optional parameter)
  //   *  \param Strategy* : the Strategy linked to the Model (optional parameter)
  //   */
  //  Model(char *xmlFile, float t=-1, float t0=-1, float T=-1, NonSmoothDynamicalSystem* nsds=NULL, Strategy* strategy=NULL);

  /** \fn Model(char *xmlFile, float t0, float T, string title, string author, string description, string date, string xmlSchema)
   *  \brief allows to create the Model with an xml file, or the needed data
   *  \param char * : the input XML file (optional parameter)
   *  \param float : the value for t0 (optional parameter)
   *  \param float : the value for T (optional parameter)
   *  \param string : the title of the Model (optional parameter)
   *  \param string : the author of the Model (optional parameter)
   *  \param string : the description of the Model (optional parameter)
   *  \param string : the date of the Model (optional parameter)
   *  \param string : the xml schema of the Model (optional parameter)
   *  \exception RuntimeException
   */
  Model(char *xmlFile, /*float t = -1.0,*/ float t0 = -1.0, float T = -1.0,
        string title = "title", string author = "author", string description = "description",
        string date = "date", string xmlSchema =/*DEFAULT_XMLSCHEMA*/ "none");

  /** \fn Model(cfloat t0, float T, string title, string author, string description, string date, string xmlSchema)
   *  \brief allows to create the Model with an xml file, or the needed data
   *  \param float : the value for t0 (optional parameter)
   *  \param float : the value for T (optional parameter)
   *  \param string : the title of the Model (optional parameter)
   *  \param string : the author of the Model (optional parameter)
   *  \param string : the description of the Model (optional parameter)
   *  \param string : the date of the Model (optional parameter)
   *  \param string : the xml schema of the Model (optional parameter)
   *  \exception RuntimeException
   */
  Model(float t0, float T,
        string title = "title", string author = "author", string description = "description",
        string date = "date", string xmlSchema =/*DEFAULT_XMLSCHEMA*/ "none");

  ~Model();

  // getter/setter
  /** \fn inline string getTitle()
   *  \brief allows to get the title of the simulation
   *  \return string : the title
   */
  inline string getTitle() const
  {
    return this->title;
  }

  /** \fn inline string getAuthor()
   *  \brief allows to get the author of the simulation
   *  \return string : the author
   */
  inline string getAuthor() const
  {
    return this->author;
  }

  /** \fn inline string getDescription()
   *  \brief allows to get the description of the simulation
   *  \return string : the description
   */
  inline string getDescription() const
  {
    return this->description;
  }

  /** \fn inline string getDate()
   *  \brief allows to get the date of the simulation
   *  \return string : the date
   */
  inline string getDate() const
  {
    return this->date;
  }

  /** \fn inline string getXMLSchema()
   *  \brief allows to get the xml schema of the simulation
   *  \return string : the xml schema
   */
  inline string getXMLSchema() const
  {
    return this->xmlSchema;
  }

  /** \fn inline void setTitle(string s)
   *  \brief allows to set the title of the simulation
   *  \param string : the title
   */
  inline void setTitle(const string s)
  {
    this->title = s;
  }

  /** \fn inline void setAuthor(string s)
   *  \brief allows to set the author of the simulation
   *  \param string : the author
   */
  inline void setAuthor(const string s)
  {
    this->author = s;
  }

  /** \fn inline void setDescription(string s)
   *  \brief allows to set the description of the simulation
   *  \param string : the description
   */
  inline void setDescription(const string s)
  {
    this->description = s;
  }

  /** \fn inline void setDate(string s)
   *  \brief allows to set the date of the simulation
   *  \param string : the date
   */
  inline void setDate(const string s)
  {
    this->date = s;
  }

  /** \fn void setXMLSchema(string)
   *  \brief allows to set the XML Schema to use to parse XML file of the platform
   *  \param string : the name of the XML Schema file
   */
  void setXMLSchema(const string str);


  /** \fn inline double getCurrentT()
   *  \brief allows to get the current time of the simulation : t
   *  \return double : the value of the curent time
   */
  inline double getCurrentT() const
  {
    return this->t;
  }

  /** \fn inline double getT0()
   *  \brief allows to get the the initial time  of the simulation: t0
   *  \return double : the value of the initial time
   */
  inline double getT0() const
  {
    return this->t0;
  }

  /** \fn inline double getFinalT()
   *  \brief allows to get the final time of the simulation : T
   *  \return double : the value of the final time
   */
  inline double getFinalT() const
  {
    return this->T;
  }


  /** \fn inline void setCurrentT(double)
   *  \brief allows to set the current time of the simulation
   *  \param double t : the value to set the current time
   */
  inline void setCurrentT(const double t)
  {
    this->t = t;
  }

  /** \fn inline void setT0(double)
   *  \brief allow to set the initial time  of the simulation
   *  \param double t0: the value to set the initial time
   */
  inline void setT0(const double t0)
  {
    this->t0 = t0;
  }

  /** \fn inline void setFinalT(double)
   *  \brief allows to set the final time
   *  \param double T : the value to set the final time
   */
  inline void setFinalT(const double T)
  {
    this->T = T;
  }

  /** \fn inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystem(void)
   *  \brief allows to get the NonSmoothDynamicalSystem of the Model
   *  \return the NonSmoothDynamicalSystem of the Model
   */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystem(void) const
  {
    return this->nsds;
  }

  /** \fn inline Strategy getStrategy(void)
   *  \brief allows to get the Strategy of the Model
   *  \return the Strategy of the Model
   */
  inline Strategy* getStrategy(void) const
  {
    return this->strategy;
  }


  /** \fn inline void setNonSmoothDynamicalSystem(NonSmoothDynamicalSystem *nsds)
   *  \brief allows to set the NonSmoothDynamicalSystem of the Model
   *  \param NonSmoothDynamicalSystem nsds : the NonSmoothDynamicalSystem to set
   */
  inline void setNonSmoothDynamicalSystem(NonSmoothDynamicalSystem *nsds)
  {
    this->nsds = nsds;
  };

  /** \fn inline void setStrategy(Strategy *str)
   *  \brief allows to set the Strategy of the Model
   *  \param Strategy *str : the Strategy to set
   */
  inline void setStrategy(Strategy *str)
  {
    this->strategy = str;
  };

  /** \fn inline SiconosModelXML getSiconosModelXML()
   *  \brief allows to get the SiconosModelXML of the Model
   *  \return SiconosModelXML : the object SiconosModelXML of the Model
   */
  inline SiconosModelXML* getSiconosModelXML() const
  {
    return this->modelxml;
  }

  /** \fn inline void setSiconosModelXML(SiconosModelXML *modelxml)
   *  \brief allows to set the SiconosModelXML of the Model
   *  \param SiconosModelXML *modelxml : the SiconosModelXML to set
   */
  inline void setSiconosModelXML(SiconosModelXML *modelxml)
  {
    this->modelxml = modelxml;
  }

  /////////////////////////////

  /** \fn bool isModelComplete(void)
   *  \brief determines if there are enough data to formalise and solve the problem
   *  \return tru if the Model is complete
   */
  bool isModelComplete(void);

  /** \fn readModel(char*)
   *  \brief reads input file to fill the formalisation of the system
   *  \param char* : the data file which must be read
   */
  void readModel(char*);

  /** \fn saveToXMLFile(char*)
   *  \brief saves into output file the data of the system
   *  \param char* : the data file which must be written
   */
  void saveToXMLFile(char*);

  /** \fn saveToDOMTree()
   *  \brief saves into the DOM tree all the data of the system
   */
  void saveToDOMTree();

  /** \fn void runSimulation(void)
   *  \brief launches the simulation
   */
  void runSimulation(void);

  /** \fn void doOneStep(void)
   *  \brief makes the simulation go on in time
   */
  void doOneStep(void);

  /** \fn void savePlatformToXML()
   *  \brief copy the data of the plateform to the XML DOM tree
   *  \exception RuntimeException
   */
  void savePlatformToXML();

  //  /** \fn void createModel(char *xmlFile, float t0, float T, string title, string author, string description, string date, string xmlSchema)
  //   *  \brief allows to create the Model with an xml file, or the needed data
  //   *  \param char * : the input XML file (optional parameter)
  //   *  \param float : the value for t0 (optional parameter)
  //   *  \param float : the value for T (optional parameter)
  //   *  \param string : the title of the Model (optional parameter)
  //   *  \param string : the author of the Model (optional parameter)
  //   *  \param string : the description of the Model (optional parameter)
  //   *  \param string : the date of the Model (optional parameter)
  //   *  \param string : the xml schema of the Model (optional parameter)
  //   *  \exception RuntimeException
  //   */
  //  void createModel(char *xmlFile=NULL, /*float t = -1.0,*/ float t0 = -1.0, float T = -1.0,
  //          string title="title", string author="author", string description="description",
  //          string date="date", string xmlSchema=/*DEFAULT_XMLSCHEMA*/ "none");
  //
  //  /** \fn void createModel(char *xmlFile, float t0, float T, string title, string author, string description, string date, string xmlSchema)
  //   *  \brief allows to create the Model with an xml file, or the needed data
  //   *  \param float : the value for t0 (optional parameter)
  //   *  \param float : the value for T (optional parameter)
  //   *  \param string : the title of the Model (optional parameter)
  //   *  \param string : the author of the Model (optional parameter)
  //   *  \param string : the description of the Model (optional parameter)
  //   *  \param string : the date of the Model (optional parameter)
  //   *  \param string : the xml schema of the Model (optional parameter)
  //   *  \exception RuntimeException
  //   */
  //  void createModel(float t0 /*= -1.0*/, float T /*= -1.0*/,
  //          string title="title", string author="author", string description="description",
  //          string date="date", string xmlSchema=/*DEFAULT_XMLSCHEMA*/ "none");

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

  /** \fn void display()
   *  \brief display the data of the Model
   */
  void display() const ;

  /** \fn int xmlSchemaValidated(string xmlFile, string xmlSchema)
   *  \brief checks if the xmlFile given respects the xmlSchema given
   *  \param string : the xml input file to check
   *  \param string : the xml schema
   *  \return int : 1 if the xml file respects the schema
   */
  int xmlSchemaValidated(string xmlFile, string xmlSchema = "");


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
  Strategy* createStrategy(string type);

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



  /******************************
   * Configuration functions :
   ******************************/

  /** \fn void setMatrixMaxSize( int max )
   *  \brief allows to change the value MatrixMaxSize which determines the bound over between
   * external file storing and xml input/outpout file storing
   *  \param int :  the value to assign to MatrixMaxSize
   */
  void setMatrixMaxSize(const int max);

  /** \fn void setVectorMaxSize( int max )
   *  \brief allows to change the value VectorMaxSize which determines the bound over between
   * external file storing and xml input/outpout file storing
   *  \param int :  the value to assign to VectorMaxSize
   */
  void setVectorMaxSize(const int max);

  /** \fn void setFileStorage( string fs )
   *  \brief allows to change the value fileStorage which determines the format
   * of external file save (binary or ascii)
   *  \param string :  the value to assign to fileStorage
   */
  void setFileStorage(const string fs);


protected:
  /** \fn linkModelXML
   *  \brief set the link between each entities of the plateform and the ...XML object corresponding
   */
  void linkModelXML(void);

  /** \fn void fillModelWithModelXML()
   *  \brief uses the ModelXML of the Model to fill the fields of the Model
   *  \exception RuntimeException
   */
  void fillModelWithModelXML();


private:
  /** contains the current time of the simulation */
  double t;

  /** contains the time of the start of the simulation */
  double t0;

  /** contains the time of the end of the simulation */
  double T;

  /** contains the strategy to solve the NonSmoothDynamicalSystem */
  Strategy *strategy;

  /** contains the NonSmoothDynamicalSystem of the simulation */
  NonSmoothDynamicalSystem * nsds;

  /** XML object linked to the Model */
  SiconosModelXML *modelxml;

  /** information about the Model */
  string title, author, description, date, xmlSchema;
};

#endif // MODEL_H

