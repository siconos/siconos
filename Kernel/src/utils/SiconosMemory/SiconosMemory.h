//$Id: SiconosMemory.h,v 1.8 2005/02/11 17:36:06 charlety Exp $

/** \class SiconosException
*   \brief This class allowa to store vectors of previous steps of the simulation
*   \author JB Charlety
*   \version 1.0
*   \date (Creation) 07/06/2004
*/

#ifndef SICONOSMEMORY_H
#define SICONOSMEMORY_H

#include <vector>
//#include "SiconosVector.h"
#include "NewSiconosVector.h"
#include "SiconosMemoryException.h"
#include "SiconosMemoryXML.h"

using namespace std;

class SiconosMemory
{

public:

  /**
   * \fn SiconosMemory()
   * \brief basic constructor.
   * memorySize and nbVectorsinMemory are set to 0. No allocation of memory for the pointers on siconosVector of vectorMemory.
   */
  SiconosMemory();

  /**
   * \fn SiconosMemory(SiconosMemoryXML *)
   * \brief constructor with the DOM tree.
   * \param SiconosMemoryXML * : the XML object which contains the data of the memory
   * memorySize is set with the parameter, and the memory is allocated for this number of SiconosVector
   */
  SiconosMemory(SiconosMemoryXML *);

  /**
   * \fn SiconosMemory(int)
   * \brief constructor with size parameter.
   * \param int : the size of the memory
   * memorySize is set with the parameter, and the memory is allocated for this number of SiconosVector
   */
  SiconosMemory(int);

  /**
   * \fn SiconosMemory(int)
   * \brief constructor with size parameter and the DOM tree.
   * \param int : the size of the memory
   * \param SiconosMemoryXML * : the XML object which contains the data of the memory
   * memorySize is set with the parameter, and the memory is allocated for this number of SiconosVector
   * and reload the data of the DOM tree.
   * This constructor uses the data of the SiconosMemoryXML but the maxSize of the SiconosMemory will be the one defined in the first parameter
   */
  SiconosMemory(int, SiconosMemoryXML *);

  /**
   * \fn SiconosMemory(vector<SiconosVector*>)
   * \brief constructor with vector parameter.
   * \param vector<SiconosVector*> : the vector of siconosVector which must be stored
   * memorySize is set to the size of the vector given in parameters
   */
  SiconosMemory(vector<SiconosVector*>);

  /**
   * \fn SiconosMemory(int, vector<SiconosVector*>)
   * \brief constructor with size and vector parameter.
   * \param int : the size of the memory
   * \param vector<SiconosVector*> : the vector of siconosVector which must be stored
   * this constructor is useful if the vector given in parameters has a size lower than the normal size of the memory
   */
  SiconosMemory(int, vector<SiconosVector*>);

  /**
   * \fn SiconosMemory(SiconosMemory&)
   * \brief constructor by copy
   */
  SiconosMemory(const SiconosMemory&);

  /**
   * \fn ~SiconosMemory()
   * \brief destructor
   * delete the SiconosVectors allocated in vectorMemory if memorySize > 0
   */
  ~SiconosMemory();

  /*************************************************************************/

  /**
   * \fn int getMemorySize()
   * \brief gives the size of the memory
   * \return int >= 0
   */
  inline int getMemorySize() const
  {
    return memorySize;
  };

  /**
   * \fn void setSiconosMemorySize( int max )
   * \brief set the max size of the SiconosMemory
   * \param int : the max size for this SiconosMemory
   */
  inline void setSiconosMemorySize(int max)
  {
    this->memorySize = max;
  };

  /**
   * \fn int getNbVectorsInMemory()
   * \brief gives the numbers of SiconosVectors currently stored in the memory
   * \return int >= 0
   */
  inline int getNbVectorsInMemory() const
  {
    return nbVectorsInMemory;
  };

  /**
   * \fn vector<SiconosVector*> getVectorMemory()
   * \brief gives the vector of SiconosVectors of the memory
   * \return stl vector od siconosVector
   */
  inline vector<SiconosVector*> getVectorMemory() const
  {
    return vectorMemory;
  };

  /**
   * \fn SiconosVector* getSiconosVector(int)
   * \brief gives a SiconosVectors of the memory
   * \param int the position in the memory of the wanted SiconosVector
   * \return SiconosVector* if the parameter has its value in [0, nbVectorsInMemory[
   */
  SiconosVector* getSiconosVector(int) const;

  /**
   * \fn void setVectorMemory(vector<SiconosVector*>)
   * \brief fill the memory with a vector of siconosVector
   * \param vector<SiconosVector*>
   * memorySize is set to the size of the vector given in parameters
   */
  void setVectorMemory(vector<SiconosVector*>);

  /**
   * \fn void setVectorMemory(vector<SiconosVector*>)
   * \brief fill the memory with a vector of siconosVector
   * \pram int : the size of the memory
   * \param vector<SiconosVector*>
   * this function is useful if the vector given in parameters has a size lower than the normal size of the memory
   */
  void setVectorMemory(int, vector<SiconosVector*>);

  /** \fn inline SiconosModelXML getSiconosModelXML()
   *  \brief allows to get the SiconosMemoryXML of the SiconosMemory
   *  \return SiconosMemoryXML* : the object SiconosMemoryXML of the SiconosMemory
   */
  inline SiconosMemoryXML* getSiconosMemoryXML()
  {
    return this->memoryXML;
  }

  /**
   * \fn void swap(SiconosVector*)
   * \brief puts a SiconosVector in the memory
   * \param SiconosVector* : the SiconosVector we want to put in memory
   * this function moves older SiconosVectors of vectorMemory of one position (this operation handles only pointers)
   * and COPIES (the values of) the SiconosVector given in parameter in the position 0 of vectorMemory
   */
  void swap(SiconosVector*);

  /**
   * \fn void swap(SiconosVector*)
   * \brief displays the data of the memory object
   */
  void display() const;

  /** \fn void fillMemoryWithMemoryXML()
   *  \brief uses the SiconosMemoryXML of the SiconosMemory to fill the fields of the SiconosMemory
   *  \exception SiconosMemoryException
   */
  void fillMemoryWithMemoryXML();

  /** \fn void saveMemorySizeToXML()
   *  \brief copy the max size of the SiconosMemory to the XML DOM tree
   *  \exception SiconosMemoryException
   */
  inline void saveMemorySizeToXML()
  {
    IN("void SiconosMemory::saveMemorySizeToXML()\n");
    if (this->memoryXML != NULL) this->memoryXML->setSiconosMemorySize(this->memorySize);
    else SiconosMemoryException::selfThrow("SiconosMemory::saveMemorySizeToXML() - memoryXML object == NULL");
    OUT("void SiconosMemory::saveMemorySizeToXML()\n");
  }

  SiconosMemory& operator = (const SiconosMemory&);


private:
  /** the maximum size of the memory (i.e the max numbers of SiconosVectors it can store) */
  int memorySize;

  /** the real number of SiconosVectors currently stored in the Memory. This number is less than or equal to memorysize */
  int nbVectorsInMemory;

  /** the stl vector which contains the SiconosVectors kept in memory */
  vector<SiconosVector*> vectorMemory;

  /** link to the XML for SiconosMemory objects */
  SiconosMemoryXML * memoryXML;

};

#endif // SICONOSMEMORY_H

//$Log: SiconosMemory.h,v $
//Revision 1.8  2005/02/11 17:36:06  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.7  2004/09/10 11:26:25  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.6  2004/08/04 11:03:23  jbarbier
//- about the SiconosMemory : when a SiconosMemory has a maxSize greater than the
//number of steps in memory required by an integrator, the oldest SiconosVector
//are deleted
//
//- the way to initialize the SiconosMemory by the integrator has been updated to
//match with these changes
//
//Revision 1.5  2004/07/16 06:23:31  charlety
//
//_ modification of the operator = of SiconosVector : the left operand can now be a composite vector
//
//_ added functions to load / save a SiconosMemory in an XML file
//
//Revision 1.4  2004/07/13 12:44:53  jbarbier
//- integration of the SiconosMemory to the plateform
//- files SiconosMemoryXML added
//
//Revision 1.3  2004/07/09 11:14:53  charlety
//
//_ Added a constructor by copy and an operator = in class SiconosMemory
//_ the getters on memory in DynamicalSystems return now some pointers
//
//Revision 1.2  2004/07/08 09:13:43  charlety
//
//_ creation of a Siconos exception dedicated to the class SiconosMemory
//_ new tests for SiconosMemory
//
//Revision 1.1  2004/07/07 13:53:13  charlety
//
//_ First version of the memory object
//_ some little bugs corrected otherwhere
//