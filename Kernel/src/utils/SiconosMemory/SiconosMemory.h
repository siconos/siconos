
/** \class SiconosException
*   \brief This class allowa to store vectors of previous steps of the simulation
*  \author SICONOS Development Team - copyright INRIA
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

//using namespace std;

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

