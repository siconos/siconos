//$Id: NewSiconosVector.cpp,v 1.11 2005/02/11 13:30:40 charlety Exp $
#include "NewSiconosVector.h"

SiconosVector::SiconosVector()
{
  //IN("SiconosVector::SiconosVector() \n");

  //OUT("SiconosVector::SiconosVector() \n");
}

SiconosVector::~SiconosVector()
{
  //IN("SiconosVector::~SiconosVector() \n");

  //OUT("SiconosVector::~SiconosVector() \n");
}


SiconosVector& SiconosVector::operator = (const SiconosVector& v)
{
  //IN(" SiconosVector::operator =\n");
  cout << " SiconosVector::operator =\n";

  if (this == &v)
  {
    //OUT(" SiconosVector::operator =\n");
    return *this;
  }
  else if (this->size() == v.size())
  {
    cout << "WARNING : operator = used with different types of vectors." << endl;
    for (int i = 0; i < this->size() ; i++)
      (*this)(i) = v(i);
    return *this;
  }
  else
    SiconosVectorException::selfThrow(" SiconosVector::operator = --- case not unforseen");
  //OUT(" SiconosVector::operator =\n");
}




void SiconosVector::zero()
{
  //IN("SiconosVector::zero() \n");

  const int size = this->size();
  for (int i = 0; i < size; i++)
    (*this)(i) = 0.0;

  //OUT("SiconosVector::zero() \n");
}


string SiconosVector::toString()
{
  //IN("SiconosVector::toString() \n");

  char element[100];
  int i = 0, end = 0;
  string vectorContent = "";
  while (i < this->size())
  {
    strcpy(element, "");
    sprintf(element, N_DOUBLE_PRECISION, (*this)(i));
    if (i > 0)
    {
      vectorContent += " ";
      vectorContent += element;
    }
    else vectorContent = element;
    i++;
  }

  //OUT("SiconosVector::toString() \n");
  return vectorContent;
}




//$Log: NewSiconosVector.cpp,v $
//Revision 1.11  2005/02/11 13:30:40  charlety
//_ added or modified some doxygen comments in SiconosMatrix.
//_ the function "getArray" is vector is now pure virtual, and implemented in CompositeVector (it returns NULL if the vector is composite).
//
//Revision 1.10  2005/02/02 15:54:51  jbarbier
//- sample RollingBalls added
//
//- function getArray() added to SimpleVector to return the pointer on the array of double values
//
//Revision 1.9  2004/09/14 13:24:53  charlety
//
//_ changes in the interface of SiconosVector
//
//Revision 1.8  2004/09/09 14:32:46  charlety
//
//_ New tests for operators of multiplication between vectors and matrices.
//
//Revision 1.7  2004/08/19 15:21:27  charlety
//
//_ SimpleVector and CompositeVector in progress.
//_ for the operators, we prefer now using directly functions of Blas1++ instead
//  of these of Blas++.h
//
//Revision 1.6  2004/08/18 14:53:23  charlety
//
//_ use of Lapack routines for operations on SimpleVector (in progress)
//
//Revision 1.5  2004/08/16 09:40:01  charlety
//
//_ new tests for the simpleVector
//
//Revision 1.4  2004/08/12 11:15:30  charlety
//
//_ SimpleVector developed at 90% (the functions of computations with matrices remain to do.)
//
//Revision 1.3  2004/08/11 14:16:07  charlety
//
//_ NewSiconosVector in progress...(NewSiconosVector is an abstract class and
//  SimpleVector inherits of NewSiconosVector).
//
//Revision 1.2  2004/07/30 14:21:54  charlety
//
//_ new functions and tests for the new SiconosVector
//
//Revision 1.1  2004/07/30 12:11:15  charlety
//
//_ try of a new SiconosVector : the composite vector is now a separated class. the development is in progress. For the moment, the platform uses the old SiconosVector.
//