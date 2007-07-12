#include "SiconosMatrix.h"

// Constructor with the type-number
SiconosMatrix::SiconosMatrix(unsigned int newNum): dimRow(0), dimCol(0), num(newNum)
{}

// Constructor with the dimensions and the type-number
SiconosMatrix::SiconosMatrix(unsigned int newNum, unsigned int row, unsigned int col): dimRow(row), dimCol(col), num(newNum)
{}

BlockIterator1 SiconosMatrix::begin()
{
  SiconosMatrixException::selfThrow("SiconosMatrix::begin(): reserved to BlockMatrix");
  BlockIterator1 it;
  return it;
};

BlockIterator1 SiconosMatrix::end()
{
  SiconosMatrixException::selfThrow("SiconosMatrix::end(): reserved to BlockMatrix");
  BlockIterator1 it;
  return it;
};

ConstBlockIterator1 SiconosMatrix::begin() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::begin(): reserved to BlockMatrix");
  ConstBlockIterator1 it;
  return it;
};

ConstBlockIterator1 SiconosMatrix::end() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::end(): reserved to BlockMatrix");
  ConstBlockIterator1 it;
  return it;
};

const Index * SiconosMatrix::getTabRowPtr() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::getTabRowPtr() : not implemented for this type of matrix (Simple?) reserved to BlockMatrix.");
  // fake to avoid error on warning.
  Index * tmp = NULL;
  return tmp;
}

const Index * SiconosMatrix::getTabColPtr() const
{
  SiconosMatrixException::selfThrow("SiconosMatrix::getTabColPtr() : not implemented for this type of matrix (Simple?) reserved to BlockMatrix.");
  // fake to avoid error on warning.
  Index * tmp = NULL;
  return tmp;
}

//=====================
// matrices comparison
//=====================
const bool isComparableTo(const  SiconosMatrix* m1, const  SiconosMatrix* m2)
{
  // return:
  // - true if one of the matrices is a Simple and if they have the same dimensions.
  // - true if both are block but with blocks which are facing each other of the same size.
  // - false in other cases

  if ((!m1->isBlock() || !m2->isBlock()) && (m1->size(0) == m2->size(0)) && (m1->size(1) == m2->size(1)))
    return true;

  const Index * I1R = m1->getTabRowPtr();
  const Index * I2R = m2->getTabRowPtr();
  const Index * I1C = m1->getTabColPtr();
  const Index * I2C = m2->getTabColPtr();

  return ((*I1R == *I2R) && (*I1C == *I2C));
}

void scal(double a, const SiconosMatrix& A, SiconosMatrix& B)
{
  // To compute B = a * A

  if (&A == &B)
    B *= a;
  else
  {
    unsigned int numA = A.getNum();
    unsigned int numB = B.getNum();

    if (numB == 6 || numB == 7) // B = 0 or identity.
      SiconosMatrixException::selfThrow("scal(a,A,B) : forbidden for B being a zero or identity matrix.");

    if (numA == 6)
      B.zero();
    else if (numA == 7)
    {
      B.eye();
      B *= a;
    }
    else
    {
      if (numA == numB) // if A and B are of the same type ...
      {
        switch (numA)
        {

        case 0: // A and B are block
          if (isComparableTo(&A, &B))
          {
            BlockIterator1 itB1;
            BlockIterator2 itB2;
            ConstBlockIterator1 itA1 = A.begin();
            ConstBlockIterator2 itA2;
            for (itB1 = B.begin(); itB1 != B.end(); ++itB1)
            {
              itA2 = itA1.begin();
              for (itB2 = itB1.begin(); itB2 != itB1.end(); ++itB2)
              {
                scal(a, **itA2++, **itB2);
              }
              itA1++;
            }
          }
          else // if A and B are not "block-consistent"
          {
            for (unsigned int i = 0; i < A.size(0); ++i)
              for (unsigned int j = 0; j < A.size(1); ++j)
                B(i, j) = a * A(i, j);

          }
          break;

        case 1: // if both are dense
          noalias(*B.getDensePtr()) = a ** A.getDensePtr();
          break;
        case 2:
          noalias(*B.getTriangPtr()) = a ** A.getTriangPtr();
          break;
        case 3:
          noalias(*B.getSymPtr()) = a ** A.getSymPtr();
          break;
        case 4:
          noalias(*B.getSparsePtr()) = a ** A.getSparsePtr();
          break;
        case 5:
          noalias(*B.getBandedPtr()) = a ** A.getBandedPtr();
          break;
        }
      }
      else // if A and B are of different types.
      {
        if (numA == 0 || numB == 0) // if A or B is block
        {
          B = A;
          B *= a;
        }
        else
        {
          if (numB != 1)
            SiconosMatrixException::selfThrow("scal(a,A,B) failed. A and B types do not fit together.");

          switch (numB)
          {
          case 1:
            noalias(*B.getDensePtr()) = a ** A.getDensePtr();
            break;
          case 2:
            noalias(*B.getDensePtr()) = a ** A.getTriangPtr();
            break;
          case 3:
            noalias(*B.getDensePtr()) = a ** A.getSymPtr();
            break;
          case 4:
            noalias(*B.getDensePtr()) = a ** A.getSparsePtr();
            break;
          case 5:
            noalias(*B.getDensePtr()) = a ** A.getBandedPtr();
            break;
          }
        }
      }
    }
  }
}
