#include "SimpleVector.h"

SimpleVector::SimpleVector()
{
  IN("SimpleVector() \n");

  this->composite = false;
  lavd = LaVectorDouble();

  OUT("SimpleVector() \n");
}


SimpleVector::SimpleVector(const string file, const bool ascii)
{
  IN("SimpleVector(const string file, const bool ascii) \n");
  string mode;
  this->composite = false;
  if (ascii) mode = "ascii";
  else mode = "binary";
  this->read(file, mode);

  OUT("SimpleVector(const string file, const bool ascii) \n");
}


SimpleVector::SimpleVector(const vector<double> v)
{
  IN("SimpleVector(const vector<double> v) \n");

  this->composite = false;
  this->setValues(v);

  OUT("SimpleVector(const vector<double> v) \n");
}


SimpleVector::SimpleVector(const SiconosVector& v)
{
  IN("SimpleVector(const SiconosVector& v) \n");

  this->composite = false;
  *this = v;

  OUT("SimpleVector(const SiconosVector& v) \n");
}


SimpleVector::SimpleVector(const SimpleVector& v)
{
  IN("SimpleVector(const SiconosVector& v) \n");

  this->composite = false;
  *this = v;

  OUT("SimpleVector(const SiconosVector& v) \n");
}

SimpleVector::SimpleVector(const int size)
{
  IN("SimpleVector (const int size) \n");
  if (size < 0)
    SiconosVectorException::selfThrow(" SimpleVector::SimpleVector(const int size)   -- can't initialize a simpleVector with a negative size");

  this->composite = false;
  LaVectorDouble lvd(size);
  lavd.resize(lvd);
  lavd = lvd;

  OUT("SimpleVector (const int size) \n");
}

SimpleVector::~SimpleVector()
{
  IN("~SimpleVector() \n");

  OUT("~SimpleVector() \n");
}

/***********************************************************************************************/

void SimpleVector::display() const
{
  IN("SimpleVector::display() \n");

  cout << "| size : " << this->size() << endl;
  cout << "| isComposite : " << this->isComposite() << endl;
  if (this->size() <= M_MAXSIZEFORDISPLAY)
  {
    for (int i = 0; i < this->size(); i++)
      cout << lavd(i) << " ";
    cout << endl;
  }
  else cout << "Display SiconosVector : vector too large" << endl;

  OUT("SimpleVector::display() \n");
}


double& SimpleVector::operator()(const int unsigned index)
{
  IN("SimpleVector::operator()(int unsigned index) \n");

  if (index >= this->lavd.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator()   -- index greater than vector size");

  OUT("SimpleVector::operator()(int unsigned index) \n");
  return (this->lavd)(index);
}


double SimpleVector::operator()(const int unsigned index) const
{
  IN("SimpleVector::operator() \n");

  if (index >= this->lavd.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator()   -- index greater than vector size");

  OUT("SimpleVector::operator()\n");
  return  this->lavd(index);
}


void SimpleVector::setValues(const vector<double> v)
{
  IN(" void SimpleVector::setValues(const vector<double> v)  \n");

  if (v.size() < 0)
    SiconosVectorException::selfThrow("SimpleVector::setValues(const vector<double> v)  -- cannot set values : the size of v is negative");

  int i;
  LaVectorDouble lvd(v.size());
  for (i = 0; i < v.size(); i++)
    lvd(i) = v[i];
  this->lavd.resize(lvd);
  this->lavd = lvd;

  OUT(" void SimpleVector::setValues(const vector<double> v)  \n");
}

void SimpleVector::zero()
{
  double *array = this->lavd.addr();
  const int size = this->size();
  for (int i = 0; i < size; i++)
    array[i] = 0.0;
}


int SimpleVector::size() const
{
  IN("int SimpleVector::size() \n");
  OUT("int SimpleVector::size() \n");
  return lavd.size();
}


bool SimpleVector::read(string fileName, string mode)
{
  IN(" SimpleVector::read(string fileName, string mode) \n");

  bool res = false;
  //this->lavd = LaVectorDouble();

  if (mode == "binary")
  {
    FILE * inFile = fopen(fileName.c_str(), "rb");    // open the input file in binary mode
    if (inFile == NULL)
    {
      SiconosVectorException::selfThrow(" SimpleVector::read : Fail to open file \"" + fileName + "\"");
    }
    int size, nbVectors;
    fread((char *) &nbVectors, sizeof(int), 1, inFile);   // read nbVectors
    if (nbVectors == 1)
    {
      fread((char *) &size, sizeof(int), 1, inFile);   // read size
      if (size > 0)
      {
        LaVectorDouble lvd(size);
        for (int i = 0; i < size; i++)
          fread((char*) & ((lvd))(i), sizeof(double), 1, inFile); // read a double
        this->lavd.resize(lvd);
        this->lavd = lvd;
      }
      else if (size == 0)
      {
        this->lavd =  LaVectorDouble();
      }
      else
      {
        SiconosVectorException::selfThrow(" SimpleVector::read : try to read a vector with a negative size");
      }
    }
    else if (nbVectors > 1)
    {
      SiconosVectorException::selfThrow(" SimpleVector::read : try to read a file saved by a composite vector");
    }
    else SiconosVectorException::selfThrow(" SimpleVector::read : nbVectors < 0");

    fclose(inFile);
    res = true;
  }
  else if (mode == "ascii")
  {
    ifstream inFile(fileName.c_str(),  ifstream::in);

    if (inFile == NULL)
    {
      SiconosVectorException::selfThrow(" SimpleVector::read : Fail to open file \"" + fileName + "\"");
    }

    int size, nbVectors;
    double tmp;

    inFile >> nbVectors;
    if (nbVectors < 0)
      SiconosVectorException::selfThrow(" SimpleVector::read : try to read a vector with a negative size");

    // get values
    vector<double> vect;
    for (int cpt = 0; cpt < nbVectors; cpt ++)
    {
      inFile >> size; // read size
      if (size > 0)
      {
        LaVectorDouble lvd(size);

        for (int i = 0; i < size; i++)
        {
          inFile >> tmp; // read a double
          vect.push_back(tmp);
        }
      }
    }

    this->setValues(vect);

    inFile.close();
    res = true;
  }
  OUT("SimpleVector::read(string fileName, string mode) \n");
  return res;
}


bool SimpleVector::write(string fileName, string mode) const
{
  IN("SimpleVector::write(string fileName, string mode) \n");

  bool res = false;
  if ((mode != "binary") && (mode != "ascii"))
    SiconosVectorException::selfThrow("SimpleVector::write : unknown mode");

  // open the file
  ofstream outFile(fileName.c_str());           // don't forget to check that it opened

  if (!outFile.is_open())
    SiconosVectorException::selfThrow("SimpleVector::write : : Fail to open file \"" + fileName + "\"");

  int nbVectors = 1;

  if (mode == "binary")
  {
    outFile.write((char*)&nbVectors, sizeof(int));
    int size = this->size();
    outFile.write((char*)&size, sizeof(int));
    for (int i = 0; i < size; i++)
      outFile.write((char*) & (this->lavd)(i), sizeof(double));
    res = true;
  }
  else if (mode == "ascii")
  {
    outFile << nbVectors;

    int size = this->size();
    outFile << endl << size << endl;
    for (int i = 0; i < size; i++)
    {
      /* WARNING this buffer is dangerous, the size is machine dependant*/
      char buffer[30];
      sprintf(buffer, "%1.17e ", (*this)(i));
      outFile << buffer;
    }
    res = true;
  }
  outFile.close();

  OUT("SimpleVector::write(string fileName, string mode) \n");
  return res;
}

double* SimpleVector::getArray()
{
  return this->lavd.addr();
}

/*******************************************************************************
*         GENERIC INTERNAL OPERATORS                                 *
*******************************************************************************/

SimpleVector &SimpleVector::operator+=(const SiconosVector &v)
{
  IN(" SimpleVector::operator+=(const SiconosVector &) \n");
  cout << "SimpleVector::operator+=(const SiconosVector &)" << endl;
  if (this->size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator+=   -- the vectors have not the same size");

  /*  Basic implementation*/

  *this = /* *this + v*/ this->addition(v);

  OUT(" SimpleVector::operator+=(const SiconosVector &) \n");
  return *this;
}


SimpleVector &SimpleVector::operator-=(const SiconosVector &v)
{
  IN(" SimpleVector::operator-=(const SiconosVector &) \n");
  //cout<<"SimpleVector::operator-=(const SiconosVector &)"<<endl;

  if (this->size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator-=   -- the vectors have not the same size");

  /*  Basic implementation*/
  *this = /* *this - v */ this->subtraction(v);

  OUT(" SimpleVector::operator-=(const SiconosVector &) \n");
  return *this;
}


SimpleVector& SimpleVector::operator = (const SiconosVector& v)
{
  IN("SimpleVector::operator = GENERIC\n");

  //cout<<"SimpleVector::operator = GENERIC\n";
  if (this != &v)
  {
    // avoid the creation of a new LaVectorDouble
    const SimpleVector *simple = dynamic_cast<const SimpleVector*>(&v);
    if (simple != NULL)
    {
      *this = *simple; // use the specific operator
    }
    else
    {
      const int vSize = v.size();
      LaVectorDouble lvd(vSize);
      lavd.resize(lvd);
      lavd = lvd;
      for (int i = 0; i < vSize; i++)
        (*this)(i) = v(i);
    }
  }

  OUT("SimpleVector::operator = GENERIC\n");
  return *this;
}


bool SimpleVector::operator == (const SiconosVector& v) const
{
  IN(" SimpleVector::operator == \n");

  bool res = true;
  if (this != &v)
  {
    const int size = this->size();
    if (size == v.size())
    {
      if (size == 0)
        res = true;
      else for (int i = 0; i < size; i++)
        {
          if ((*this)(i) != v(i))
          {
            res = false;
            break;
          }
        }
    }
    else res =  false;
  }

  OUT(" SimpleVector::operator == \n");
  return res;
}



bool SimpleVector::operator != (const SiconosVector& v) const
{
  IN(" SimpleVector::operator != \n");

  return !(*this == v);

  OUT(" SimpleVector::operator != \n");
}


/*******************************************************************************
*         SPECIFIC INTERNAL OPERATORS                                *
*******************************************************************************/
SimpleVector& SimpleVector::operator = (const SimpleVector& v)
{
  /*
   * Could be faster without the constness of v
   */

  IN("SimpleVector::operator = SPECIFIC\n");
  //  cout<<"SimpleVector::operator = SPECIFIC\n";
  if (this != &v)
  {
    /* BLAS - LAPACK routine */
    //    cout<<"call of Lapack routine for affectation of simple vector"<<endl;
    this->lavd.copy(v.lavd);
  }

  OUT("SimpleVector::operator = SPECIFIC\n");
  return *this;

}


SimpleVector &SimpleVector::operator+=(const SimpleVector &v)
{
  IN(" SimpleVector::operator+=(const SimpleVector &) \n");
  //cout<<" SimpleVector::operator+=(const SimpleVector &) \n";

  if (this->size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator+=   -- the vectors have not the same size");

  /* BLAS - LAPACK routine */
  //this->lavd = this->lavd + v.lavd;
  Blas_Add_Mult(this->lavd, 1.0, v.lavd);

  OUT(" SimpleVector::operator+=(const SimpleVector &) \n");
  return *this;

}


SimpleVector &SimpleVector::operator-=(const SimpleVector &v)
{
  IN(" SimpleVector::operator-=(const SimpleVector &) \n");
  //cout<<" SimpleVector::operator-=(const SimpleVector &) \n";

  if (this->size() != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::operator-=   -- the vectors have not the same size");

  /* BLAS - LAPACK routine */

  Blas_Add_Mult(this->lavd, -1.0, v.lavd);

  OUT(" SimpleVector::operator-=(const SimpleVector &) \n");
  return *this;
}


SimpleVector &SimpleVector::operator*=(const double d)
{
  IN(" SimpleVector::operator*=(const double d) \n");

  /* BLAS - LAPACK routine */

  //  cout<<"call of Lapack routine for multiplication of simple vector"<<endl;
  this->lavd = this->lavd * d;

  OUT(" SimpleVector::operator*=(const double d) \n");
  return *this;
}


SimpleVector &SimpleVector::operator/=(const double d)
{
  IN(" SimpleVector::operator/=(const double d) \n");

  if (d == 0)
    SiconosVectorException::selfThrow(" SimpleVector::operator/ =   --division by 0");

  /* BLAS - LAPACK routine */

  //  cout<<"call of Lapack routine for division of simple vector"<<endl;
  this->lavd = this->lavd * (1 / d);

  OUT(" SimpleVector::operator/=(const double d) \n");
  return *this;
}

/*******************************************************************************
*         GENERIC EXTERNAL OPERATORS                                 *
/******************************************************************************/

SimpleVector operator + (const SiconosVector& v1, const SiconosVector& v2)
{
  IN(" SimpleVector operator +  (const SimpleVector& , const SiconosVector&) GENERIC \n");

  const int svSize = v1.size();

  if (svSize != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator + --- the vectors have not the same size");

  /*  Basic implementation*/
  cout << "WARNING - generic addition of SiconosVector" << endl;

  SimpleVector sv(svSize);

  //sv = v1.addition(v2);
  for (int i = 0; i < svSize; i++) sv(i) = v1(i) + v2(i);

  OUT(" SimpleVector operator +  (const SimpleVector& , const SiconosVector&) GENERIC \n");
  return sv;
}


SimpleVector operator - (const SiconosVector& v1, const SiconosVector& v2)
{
  IN(" SimpleVector operator -  (const SimpleVector& , const SiconosVector&) GENERIC \n");

  const int svSize = v1.size();

  if (svSize != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator - --- the vectors have not the same size");

  /*  Basic implementation*/
  cout << "WARNING -  generic subtraction of SiconosVector" << endl;

  SimpleVector sv(svSize);

  //sv = v1.subtraction(v2);
  for (int i = 0; i < svSize; i++) sv(i) = v1(i) - v2(i);

  OUT(" SimpleVector operator -  (const SimpleVector& , const SiconosVector&) GENERIC \n");
  return sv;
}

/*******************************************************************************
*         GENERIC INTERNAL FUNCTIONS FOR MIXED OPERATIONS            *
/******************************************************************************/

SimpleVector SimpleVector::addition(const SiconosVector& v) const
{
  IN(" SimpleVector::addition(const SiconosVector& v) GENERIC \n");

  const int svSize = this->size();

  if (svSize != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::addition(const SiconosVector& v) --- the vectors have not the same size");

  /*  Basic implementation*/
  cout << "WARNING - addition of a simpleVector with a SiconosVector" << endl;

  SimpleVector sv(svSize);
  for (int i = 0; i < svSize; i++)
    sv.lavd(i) = this->lavd(i) + v(i);

  OUT(" SimpleVector::addition(const SiconosVector& v) GENERIC \n");
  return sv;
}


SimpleVector SimpleVector::subtraction(const SiconosVector& v) const
{
  IN(" SimpleVector::subtraction(const SiconosVector& v) GENERIC \n");

  const int svSize = this->size();

  if (svSize != v.size())
    SiconosVectorException::selfThrow(" SimpleVector::subtraction(const SiconosVector& v) --- the vectors have not the same size");

  /*  Basic implementation*/
  cout << "WARNING - subtraction of a simpleVector with a SiconosVector" << endl;

  SimpleVector sv(svSize);
  for (int i = 0; i < svSize; i++)
    sv.lavd(i) = this->lavd(i) - v(i);

  OUT(" SimpleVector::subtraction(const SiconosVector& v) GENERIC \n");
  return sv;
}


/*******************************************************************************
*         SPECIFIC EXTERNAL OPERATORS                                *
*******************************************************************************/
SimpleVector operator * (const SimpleVector& v, const double d)
{
  IN(" SimpleVector operator * (const SimpleVector& v, const double d)  \n");

  /**
   * WARNING : try of different implemetations. Direct call to Blas is faster
   */

  //  // uses LaVectorDouble operator * based on a simple loop
  //  SimpleVector sv(v.size());
  //  sv.lavd = v.lavd * d;
  //  return sv;

  // call of Blas_Scale : needs a copy before multiplication
  SimpleVector sv(v);
  Blas_Scale(d, sv.lavd);
  OUT(" SimpleVector operator * (const SimpleVector& v, const double d)  \n");
  return sv;

  //  // uses LaVectorDouble operator * based on a simple loop
  //  const int svsize = v.size();
  //  SimpleVector sv(svsize);
  //  for (int i = 0; i < svsize; i++) sv(i) = v(i) * d;
  //  return sv;
}


SimpleVector operator * (const double d, const SimpleVector& v)
{
  IN(" SimpleVector::operator *  (const SimpleVector& v, const double d) \n");

  SimpleVector sv(v);
  Blas_Scale(d, sv.lavd);

  OUT(" SimpleVector::operator * (const SimpleVector& v, const double d)  \n");
  return sv;
}


SimpleVector operator / (const SimpleVector& v, const double d)
{
  IN(" SimpleVector operator /  (const SimpleVector& v, const double d) \n");

  if (d == 0.0)
    SiconosVectorException::selfThrow(" SimpleVector operator / --- division by 0");

  SimpleVector sv(v);
  Blas_Scale(1 / d, sv.lavd);

  OUT(" SimpleVector operator / (const SimpleVector& v, const double d)  \n");
  return sv;
}


SimpleVector operator + (const SimpleVector& v1, const SimpleVector& v2)
{
  IN(" SimpleVector operator +  (const SimpleVector& , const SimpleVector&) \n");

  const int svSize = v1.size();

  if (svSize != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator + --- the 2 vectors have not the same size");

  SimpleVector sv(v1);
  Blas_Add_Mult(sv.lavd, 1.0, v2.lavd);

  return sv;

  OUT("SimpleVector operator +  (const SimpleVector& , const SimpleVector&)  \n");
}


SimpleVector operator - (const SimpleVector& v1, const SimpleVector& v2)
{
  IN(" SimpleVector operator -  (const SimpleVector& , const SimpleVector&) \n");

  const int svSize = v1.size();

  if (svSize != v2.size())
    SiconosVectorException::selfThrow(" SimpleVector operator - --- the 2 vectors have not the same size");

  SimpleVector sv(v1);
  Blas_Add_Mult(sv.lavd, -1.0, v2.lavd);

  OUT("SimpleVector operator -  (const SimpleVector& , const SimpleVector&)  \n");
  return sv;
}


SimpleVector operator * (/*const*/ SiconosMatrix &m, /*const*/ SimpleVector &v)
{
  SimpleVector sv(m.size(0));
  Blas_Mat_Vec_Mult(m.getLaGenMatDouble(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}


SimpleVector matTransVecMult(SiconosMatrix &m, SimpleVector &v)
{
  SimpleVector sv(m.size(1));
  Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), v.lavd, sv.lavd, 1.0, 0.0);
  return sv;
}



