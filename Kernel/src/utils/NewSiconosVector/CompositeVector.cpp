
#include "CompositeVector.h"

CompositeVector::CompositeVector()
{
  //IN("CompositeVector() \n");

  this->composite = true;
  svref.clear();
  tabindex.clear();

  //OUT("CompositeVector() \n");
}

CompositeVector::CompositeVector(const string file, const bool ascii)
{
  //IN("CompositeVector(const string file, const bool ascii) \n");
  string mode;
  this->composite = true;
  if (ascii) mode = "ascii";
  else mode = "binary";
  this->read(file, mode);

  //OUT("CompositeVector(const string file, const bool ascii) \n");
}

CompositeVector::CompositeVector(const vector<double> v)
{
  //IN("CompositeVector(const vector<double> v) \n");

  SiconosVectorException::selfThrow(" CompositeVector::CompositeVector(const vector<double> v) -- operation available only for simpleVector");

  //OUT("CompositeVector(const vector<double> v) \n");
}

CompositeVector::CompositeVector(const SiconosVector& v)
{
  //IN("CompositeVector(const SiconosVector& v) \n");

  this->composite = true;
  this->add(v);

  //OUT("CompositeVector(const SiconosVector& v) \n");
}


CompositeVector::CompositeVector(const CompositeVector& v)
{
  //IN("CompositeVector(const SiconosVector& v) \n");

  this->composite = true;

  if (this != &v)
  {
    this->svref.clear();
    this->tabindex.clear();

    this->svref = v.svref;
    this->tabindex = v.tabindex;

    //*this = v;
  }
  else
    SiconosVectorException::selfThrow("CompositeVector(const SiconosVector& v)   -- recursive composition is not allowed");

  //OUT("CompositeVector(const SiconosVector& v) \n");
}


CompositeVector::CompositeVector(const int size)
{
  //IN("CompositeVector (const int size) \n");

  SiconosVectorException::selfThrow(" CompositeVector::CompositeVector(const int size)   -- can't initialize a CompositeVector with a  size");

  //OUT("CompositeVector (const int size) \n");
}


CompositeVector::~CompositeVector()
{
  //IN("~CompositeVector() \n");

  svref.clear();
  tabindex.clear();

  //OUT("~CompositeVector() \n");
}


void CompositeVector::display() const
{
  //IN("CompositeVector::display() \n");

  cout << "| size : " << this->size() << endl;
  cout << "| isComposite : " << this->isComposite() << endl;

  for (int i = 0; i < this->svref.size(); i ++)
  {
    cout << "\n| index : " << this->tabindex[i] << endl;
    this->svref[i]->display();
  }

  //OUT("CompositeVector::display() \n");
}


double& CompositeVector::operator()(const int unsigned index)
{
  //IN("CompositeVector::operator()(int unsigned index) \n");
  if (index >   this->size())
    SiconosVectorException::selfThrow(" CompositeVector::operator()   -- out of range");

  int indexVect = 0, pos = index, sizeTabIndex = this->tabindex.size();;

  while ((indexVect < sizeTabIndex) && (this->tabindex[indexVect] <= index))
  {
    pos = index - this->tabindex[indexVect] ;
    indexVect++;
  }

  //OUT("CompositeVector::operator()(int unsigned index) \n");
  return (*(this->svref[indexVect]))(pos);
}


double CompositeVector::operator()(const int unsigned index) const
{
  //IN("CompositeVector::operator() \n");
  //cout<<"CompositeVector::operator() \n";

  if (index >   this->size())
    SiconosVectorException::selfThrow(" CompositeVector::operator()   -- out of range");

  int indexVect = 0, pos = index, sizeTabIndex = this->tabindex.size();

  while ((indexVect < sizeTabIndex) && (this->tabindex[indexVect] <= index))
  {
    pos = index - this->tabindex[indexVect] ;
    indexVect++;
  }

  //OUT("CompositeVector::operator()\n");
  return (*(this->svref[indexVect]))(pos);
}


void CompositeVector::add(const SiconosVector& v)
{
  //IN("CompositeVector::add(const SiconosVector& v)  \n");

  SiconosVector *sv = const_cast<SiconosVector*>(&v);
  this->svref.push_back(sv);
  if (this->tabindex.size() > 0)
    this->tabindex.push_back(this->tabindex[this->tabindex.size() - 1] + v.size());
  else
    this->tabindex.push_back(v.size());


  //OUT("CompositeVector::add(const SiconosVector& v)  \n");
}


void CompositeVector::setValues(const vector<double> v)
{
  //IN(" void CompositeVector::setValues(const vector<double> v)  \n");

  SiconosVectorException::selfThrow("CompositeVector::setValues  -- this operation is not available for compositeVector");

  //OUT(" void CompositeVector::setValues(const vector<double> v)  \n");
}


int CompositeVector::size() const
{
  //IN("int CompositeVector::size() \n");

  int res = 0;
  if (this->tabindex.size() > 0)  res = this->tabindex[this->tabindex.size() - 1];
  return res;

  //OUT("int CompositeVector::size() \n");
}


bool CompositeVector::read(string fileName, string mode)
{
  //IN(" CompositeVector::read(string fileName, string mode) \n");

  bool res = false;
  //this->lavd = LaVectorDouble();

  if (mode == "binary")
  {
    FILE * inFile = fopen(fileName.c_str(), "rb");    // open the input file in binary mode
    if (inFile == NULL)
    {
      SiconosVectorException::selfThrow(" CompositeVector::read : Fail to open file \"" + fileName + "\"");
    }
    int size, nbVectors = 0, savedSize = 0;
    char essai[100];
    fread((char *) &nbVectors, sizeof(int), 1, inFile);   // read nbVectors
    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      fread((char *) &size, sizeof(int), 1, inFile);   // read size
      if (size != this->svref[cpt]->size())
        SiconosVectorException::selfThrow(" CompositeVector::read : sub-vector has a bad size.");
      for (int i = 0; i < size; i++)
        fread((char*) & (*this)(savedSize + i), sizeof(double), 1, inFile); // read a double
      savedSize += size;
    }

    fclose(inFile);
    res = true;
  }
  else if (mode == "ascii")
  {
    ifstream inFile(fileName.c_str(),  ifstream::in);

    if (inFile == NULL)
    {
      SiconosVectorException::selfThrow(" CompositeVector::read : Fail to open file \"" + fileName + "\"");
    }

    int size, nbVectors, savedSize = 0;
    double tmp;

    inFile >> nbVectors;

    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      inFile >> size; // read size
      if (size > 0)
      {
        if (size != this->svref[cpt]->size())
          SiconosVectorException::selfThrow(" CompositeVector::read : sub-vector has a bad size.");
        for (int i = 0; i < size; i++)
        {
          inFile >> tmp; // read a double
          (*this)(i + savedSize) = tmp;
        }
        savedSize += size;
      }
    }

    inFile.close();
    res = true;
  }

  return true;
  //OUT("CompositeVector::read(string fileName, string mode) \n");
}


bool CompositeVector::write(string fileName, string mode) const
{
  //IN("CompositeVector::write(string fileName, string mode) \n");

  bool res = false;
  if ((mode != "binary") && (mode != "ascii"))
    SiconosVectorException::selfThrow("SimpleVector::write : unknown mode");

  // open the file
  ofstream outFile(fileName.c_str());           // don't forget to check that it opened

  if (!outFile.is_open())
    SiconosVectorException::selfThrow("SimpleVector::write : : Fail to open file \"" + fileName + "\"");

  int nbVectors = this->svref.size();

  if (mode == "binary")
  {
    outFile.write((char*)&nbVectors, sizeof(int));
    int size = 0;
    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      size = this->svref[cpt]->size();
      outFile.write((char*)&size, sizeof(int));
      for (int i = 0; i < size; i++)
        outFile.write((char*) & (*(this->svref[cpt]))(i), sizeof(double));
    }
    res = true;
  }
  else if (mode == "ascii")
  {
    outFile << nbVectors << endl;

    int size = 0;
    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      size = this->svref[cpt]->size();
      outFile << size << endl;
      for (int i = 0; i < size; i++)
      {
        /* WARNING this buffer is dangerous, the size is machine dependant*/
        char buffer[30];
        sprintf(buffer, "%1.17e ", (*(this->svref[cpt]))(i));
        outFile << buffer;
      }
      outFile << endl;
    }
    res = true;
  }
  outFile.close();

  //OUT("CompositeVector::write(string fileName, string mode) \n");
  return res;
}


double* CompositeVector::getArray()
{
  cout << "The SiconosVector is Composite, this function is not available" << endl;
  return NULL;
}


/*******************************************************************************
*         GENERIC INTERNAL OPERATORS                                 *
*******************************************************************************/

CompositeVector &CompositeVector::operator+=(const SiconosVector &v)
{
  //IN(" CompositeVector::operator+=(const SiconosVector &) \n");
  //cout<<" CompositeVector::operator+=(const SiconosVector &) \n";

  const int size = this->size();
  if (size != v.size())
    SiconosVectorException::selfThrow(" CompositeVector::operator+=   -- the vectors have not the same size");

  for (int i = 0; i < size; i++)
    (*this)(i) += v(i);

  //OUT(" CompositeVector::operator+=(const SiconosVector &) \n");
  return *this;
}


CompositeVector &CompositeVector::operator-=(const SiconosVector &v)
{
  //IN(" CompositeVector::operator-=(const SiconosVector &) \n");
  //cout<<" CompositeVector::operator-=(const SiconosVector &) \n";

  const int size = this->size();
  if (size != v.size())
    SiconosVectorException::selfThrow(" CompositeVector::operator-=   -- the vectors have not the same size");

  for (int i = 0; i < size; i++)
    (*this)(i) -= v(i);

  //OUT(" CompositeVector::operator-=(const SiconosVector &) \n");
  return *this;
}


CompositeVector& CompositeVector::operator = (const SiconosVector& v)
{
  //IN("CompositeVector::operator = \n");

  if (this != &v)
  {
    const int size = this->size();
    if (size != v.size())
      SiconosVectorException::selfThrow(" CompositeVector::operator = GENERIC  -- the vectors have not the same size");
    else
    {
      for (int i = 0; i < size; i++)
        (*this)(i) = v(i);
    }
  }
  //OUT("CompositeVector::operator = \n");
  return *this;
}


bool CompositeVector::operator == (const SiconosVector& v) const
{
  //IN(" CompositeVector::operator == \n");
  //cout<<" CompositeVector::operator == \n";
  bool res = true;
  int i = 0;
  if (this != &v)
  {
    const int size = this->size();
    if (size == v.size())
    {
      if (size == 0)
        res = true;
      else for (i = 0; i < size; i++)
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

  //OUT(" CompositeVector::operator == \n");
  return res;
}


bool CompositeVector::operator != (const SiconosVector& v) const
{
  //IN(" CompositeVector::operator != \n");

  return !(*this == v);

  //OUT(" CompositeVector::operator != \n");
}


/******************************************************************************/
//          SPECIFIC INTERNAL OPERATORS
/******************************************************************************/

CompositeVector &CompositeVector::operator*=(const double d)
{
  //IN(" CompositeVector::operator*=(const double d) \n");

  *this = *this * d;

  //OUT(" CompositeVector::operator*=(const double d) \n");
  return *this;
}


CompositeVector &CompositeVector::operator/=(const double d)
{
  //IN(" CompositeVector::operator/=(const double d) \n");

  *this = *this / d;

  //OUT(" CompositeVector::operator/=(const double d) \n");
  return *this;
}


CompositeVector &CompositeVector::operator+=(const CompositeVector& v)
{
  //IN(" CompositeVector::operator+=(const CompositeVector& v) \n");

  *this = *this + v;

  //OUT(" CompositeVector::operator+=(const CompositeVector& v) \n");
  return *this;
}


CompositeVector &CompositeVector::operator-=(const CompositeVector& v)
{
  //IN(" CompositeVector::operator-=(const CompositeVector& v) \n");

  *this = *this - v;
  return *this;

  //OUT(" CompositeVector::operator-=(const CompositeVector& v) \n");
}


CompositeVector& CompositeVector::operator = (const CompositeVector& v)
{
  //IN("CompositeVector::operator = \n");
  //cout<<"CompositeVector::operator = SPC\n";

  if (this != &v)
  {
    {
      const int size = this->size();
      if (size != v.size())
        SiconosVectorException::selfThrow(" CompositeVector::operator = GENERIC  -- the vectors have not the same size");
      else
      {
        for (int i = 0; i < size; i++)
          (*this)(i) = v(i);
      }
    }
  }

  //OUT("CompositeVector::operator = \n");
  return *this;
}


/*******************************************************************************
*         GENERIC EXTERNAL OPERATORS                                 *
/******************************************************************************/


/*******************************************************************************
*         GENERIC EXTERNAL OPERATORS                                 *
/******************************************************************************/

CompositeVector CompositeVector::addition(const SiconosVector& v1) const
{
  //IN("CompositeVector::addition(const CompositeVector& v1)\n");

  const int size = this->size();
  const int nbVec1 = this->svref.size();
  int i = 0;

  if (size != v1.size())
    SiconosVectorException::selfThrow(" CompositeVector::addition(const CompositeVector& v1)  --- the vectors have not the same size");

  CompositeVector cv(*this);

  // the vectors have not the same structure
  cout << "WARNING : CompositeVector::addition(const CompositeVector& v1)  --- the vectors have not the same type (loops)" << endl;
  for (i = 0; i < size; i++)
    cv(i) += v1(i);

  //OUT("CompositeVector::addition(const CompositeVector& v1)\n");
  return cv;
}


CompositeVector CompositeVector::subtraction(const SiconosVector& v1) const
{
  //IN("CompositeVector::subtraction(const CompositeVector& v1)\n");

  const int size = this->size();
  const int nbVec1 = this->svref.size();
  int i = 0;

  if (size != v1.size())
    SiconosVectorException::selfThrow(" CompositeVector::subtraction(const CompositeVector& v1)  --- the vectors have not the same size");

  CompositeVector cv(*this);

  // the vectors have not the same structure
  cout << "WARNING : CompositeVector::subtraction(const CompositeVector& v1)  --- the vectors have not the same type (loops)" << endl;
  for (i = 0; i < size; i++)
    cv(i) -= v1(i);

  //OUT("CompositeVector::subtraction(const CompositeVector& v1)\n");
  return cv;
}

/*******************************************************************************
*         SPECIFIC EXTERNAL OPERATORS                                 *
/******************************************************************************/

CompositeVector operator * (const CompositeVector& v, const double d)
{
  //IN(" CompositeVector operator * (const CompositeVector& v, const double d) \n");

  const int size = v.svref.size();
  CompositeVector cv(v);
  for (int i = 0; i < size; i++)
    *(v.svref[i]) *= d;

  //OUT(" CompositeVector operator * (const CompositeVector& v, const double d)\n");
  return cv;
}


CompositeVector operator * (const double d, const CompositeVector& v)
{
  //IN(" CompositeVector friend operator * \n");

  const int size = v.svref.size();
  CompositeVector cv(v);
  for (int i = 0; i < size; i++)
    *(v.svref[i]) *= d;

  //OUT(" CompositeVector friend operator * \n");
  return cv;
}


CompositeVector operator / (const CompositeVector& v, const double d)
{
  //IN(" CompositeVector friend operator * \n");

  if (d == 0.0)
    SiconosVectorException::selfThrow(" CompositeVector operator/   --- division by 0");

  const int size = v.svref.size();
  CompositeVector cv(v);
  for (int i = 0; i < size; i++)
    *(v.svref[i]) /= d;

  //OUT(" CompositeVector friend operator * \n");
  return cv;
}


CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2)
{
  //IN("CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2) \n");

  const int size = v1.size();
  const int nbVec1 = v1.tabindex.size();
  const int nbVec2 = v2.tabindex.size();
  int i = 0;
  bool useBlas = true;

  if (size != v2.size())
    SiconosVectorException::selfThrow(" CompositeVector operator+  --- the vectors have not the same size");

  CompositeVector cv(v1);

  // comparison of tabindex
  if (nbVec1 == nbVec2)
  {
    for (i = 0; i < nbVec1; i++)
    {
      if (v1.tabindex[i] != v2.tabindex[i])
      {
        useBlas = false;
        break;
      }
    }
  }
  else useBlas = false;

  if (useBlas)
  {
    // call to BLAS function
    for (i = 0; i < nbVec1; i++)
      *(cv.svref[i]) += *(v2.svref[i]);
  }
  else
  {
    // the vectors have not the same structure
    cout << "WARNING : CompositeVector operator+  --- the vectors have not the same structure" << endl;
    for (i = 0; i < size; i++)
      cv(i) += v2(i);
  }

  //OUT("CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2) \n");
  return cv;
}


CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2)
{
  //IN("CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2) \n");

  const int size = v1.size();
  const int nbVec1 = v1.tabindex.size();
  const int nbVec2 = v2.tabindex.size();
  int i = 0;
  bool useBlas = true;

  if (size != v2.size())
    SiconosVectorException::selfThrow(" CompositeVector operator-  --- the vectors have not the same size");

  CompositeVector cv(v1);

  // comparison of tabindex
  if (nbVec1 == nbVec2)
  {
    for (i = 0; i < nbVec1; i++)
    {
      if (v1.tabindex[i] != v2.tabindex[i])
      {
        useBlas = false;
        break;
      }
    }
  }
  else useBlas = false;

  if (useBlas)
  {
    // call to BLAS function
    for (i = 0; i < nbVec1; i++)
      *(cv.svref[i]) -= *(v2.svref[i]);
  }
  else
  {
    // the vectors have not the same structure
    cout << "WARNING : CompositeVector operator-  --- the vectors have not the same structure" << endl;
    for (i = 0; i < size; i++)
      cv(i) -= v2(i);
  }

  //OUT("CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2) \n");
  return cv;
}


SimpleVector operator * (/*const*/ SiconosMatrix &m, /*const*/ CompositeVector &v)
{
  SimpleVector sv(m.size(0));

  LaVectorDouble lvd(v.size()), lvd2(sv.size());

  int size = v.size();

  for (int i = 0; i < size; i++)
  {
    // copy values in lavd
    lvd(i) = v(i);
  }

  Blas_Mat_Vec_Mult(m.getLaGenMatDouble(), lvd, lvd2, 1.0, 0.0);
  for (int i = 0; i < sv.size(); i++)
  {
    sv(i) = lvd2(i);
  }
  //cout<<"operator * (/*const*/ SiconosMatrix &m, /*const*/ SiconosVector& v) "<<endl;
  return sv;
}


SimpleVector matTransVecMult(SiconosMatrix &m, SiconosVector &v)
{
  //cout<<"SimpleVector matTransVecMult(SiconosMatrix &m, SiconosVector &v)"<<endl;

  const int size = m.size(1);

  SimpleVector sv(size);

  int i;
  LaVectorDouble lvd(v.size());
  LaVectorDouble lvd2(size);

  for (i = 0; i < v.size(); i++)
    lvd(i) = v(i);
  Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), lvd, lvd2, 1.0, 0.0);
  for (i = 0; i < size; i++)
    sv(i) = lvd2(i);

  return sv;
}


