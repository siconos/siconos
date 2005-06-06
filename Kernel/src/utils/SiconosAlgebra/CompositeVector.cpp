#include "CompositeVector.h"
using namespace std;

CompositeVector::CompositeVector()
{
  composite = true;
  svref.clear();
  tabindex.clear();
}

CompositeVector::CompositeVector(const string& file, const bool& ascii)
{
  string mode;
  composite = true;
  if (ascii) mode = "ascii";
  else mode = "binary";
  read(file, mode);
}

CompositeVector::CompositeVector(const vector<double>& v)
{
  SiconosVectorException::selfThrow(" CompositeVector::CompositeVector(const vector<double> v) -- operation available only for simpleVector");
}

CompositeVector::CompositeVector(const SiconosVector& v)
{
  composite = true;
  add(v);
}

CompositeVector::CompositeVector(const CompositeVector& v)
{
  composite = true;
  if (this != &v)
  {
    svref.clear();
    tabindex.clear();
    svref = v.svref;
    tabindex = v.tabindex;
  }
  else
    SiconosVectorException::selfThrow("CompositeVector(const SiconosVector& v)   -- recursive composition is not allowed");
}


CompositeVector::CompositeVector(const int& size)
{
  SiconosVectorException::selfThrow(" CompositeVector::CompositeVector(const int size)   -- can't initialize a CompositeVector with a  size");
}


CompositeVector::~CompositeVector()
{
  svref.clear();
  tabindex.clear();
}


void CompositeVector::display() const
{
  cout << "| size : " << size() << endl;
  cout << "| isComposite : " << isComposite() << endl;

  const int sizeV = svref.size();
  for (int i = 0; i < sizeV; i ++)
  {
    cout << "\n| index : " << tabindex[i] << endl;
    svref[i]->display();
  }
}


double& CompositeVector::operator()(const int unsigned index)
{
  //IN("CompositeVector::operator()(int unsigned index) \n");
  if (index >   size())
    SiconosVectorException::selfThrow(" CompositeVector::operator()   -- out of range");

  int indexVect = 0, pos = index, sizeTabIndex = tabindex.size();;

  while ((indexVect < sizeTabIndex) && (tabindex[indexVect] <= index))
  {
    pos = index - tabindex[indexVect] ;
    indexVect++;
  }

  //OUT("CompositeVector::operator()(int unsigned index) \n");
  return (*(svref[indexVect]))(pos);
}


double CompositeVector::operator()(const int unsigned index) const
{
  //IN("CompositeVector::operator() \n");
  //cout<<"CompositeVector::operator() \n";

  if (index >   size())
    SiconosVectorException::selfThrow(" CompositeVector::operator()   -- out of range");

  int indexVect = 0, pos = index, sizeTabIndex = tabindex.size();

  while ((indexVect < sizeTabIndex) && (tabindex[indexVect] <= index))
  {
    pos = index - tabindex[indexVect] ;
    indexVect++;
  }

  //OUT("CompositeVector::operator()\n");
  return (*(svref[indexVect]))(pos);
}


void CompositeVector::add(const SiconosVector &v)
{
  //IN("CompositeVector::add(const SiconosVector& v)  \n");

  SiconosVector *sv = const_cast<SiconosVector*>(&v);
  svref.push_back(sv);
  if (tabindex.size() > 0)
    tabindex.push_back(tabindex[tabindex.size() - 1] + v.size());
  else
    tabindex.push_back(v.size());


  //OUT("CompositeVector::add(const SiconosVector& v)  \n");
}
void CompositeVector::add(SiconosVector *v)
{
  //IN("CompositeVector::add(const SiconosVector& v)  \n");

  svref.push_back(v);
  if (tabindex.size() > 0)
    tabindex.push_back(tabindex[tabindex.size() - 1] + v->size());
  else
    tabindex.push_back(v->size());


  //OUT("CompositeVector::add(const SiconosVector& v)  \n");
}


void CompositeVector::setValues(const vector<double>& v)
{
  //IN(" void CompositeVector::setValues(const vector<double> v)  \n");

  SiconosVectorException::selfThrow("CompositeVector::setValues  -- this operation is not available for compositeVector");

  //OUT(" void CompositeVector::setValues(const vector<double> v)  \n");
}


int CompositeVector::size() const
{
  //IN("int CompositeVector::size() \n");

  int res = 0;
  if (tabindex.size() > 0)  res = tabindex[tabindex.size() - 1];
  return res;

  //OUT("int CompositeVector::size() \n");
}


bool CompositeVector::read(const string& fileName, const string& mode)
{
  //IN(" CompositeVector::read(string fileName, string mode) \n");

  bool res = false;
  //lavd = LaVectorDouble();

  if (mode == "binary")
  {
    FILE * inFile = fopen(fileName.c_str(), "rb");    // open the input file in binary mode
    if (inFile == NULL)
    {
      SiconosVectorException::selfThrow(" CompositeVector::read : Fail to open file \"" + fileName + "\"");
    }
    int sizeV, nbVectors = 0, savedSize = 0;
    char essai[100];
    fread((char *) &nbVectors, sizeof(int), 1, inFile);   // read nbVectors
    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      fread((char *) &sizeV, sizeof(int), 1, inFile);   // read size
      if (sizeV != svref[cpt]->size())
        SiconosVectorException::selfThrow(" CompositeVector::read : sub-vector has a bad size.");
      for (int i = 0; i < sizeV; i++)
        fread((char*) & (*this)(savedSize + i), sizeof(double), 1, inFile); // read a double
      savedSize += sizeV;
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

    int sizeV, nbVectors, savedSize = 0;
    double tmp;

    inFile >> nbVectors;

    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      inFile >> sizeV; // read size
      if (sizeV > 0)
      {
        if (sizeV != svref[cpt]->size())
          SiconosVectorException::selfThrow(" CompositeVector::read : sub-vector has a bad size.");
        for (int i = 0; i < sizeV; i++)
        {
          inFile >> tmp; // read a double
          (*this)(i + savedSize) = tmp;
        }
        savedSize += sizeV;
      }
    }

    inFile.close();
    res = true;
  }

  return true;
  //OUT("CompositeVector::read(string fileName, string mode) \n");
}


bool CompositeVector::write(const string& fileName, const  string& mode) const
{
  //IN("CompositeVector::write(string fileName, string mode) \n");

  bool res = false;
  if ((mode != "binary") && (mode != "ascii"))
    SiconosVectorException::selfThrow("CompositeVector::write : unknown mode");

  // open the file
  ofstream outFile(fileName.c_str());           // don't forget to check that it opened

  if (!outFile.is_open())
    SiconosVectorException::selfThrow("CompositeVector::write : : Fail to open file \"" + fileName + "\"");

  int nbVectors = svref.size();

  if (mode == "binary")
  {
    outFile.write((char*)&nbVectors, sizeof(int));
    int sizeV = 0;
    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      sizeV = svref[cpt]->size();
      outFile.write((char*)&sizeV, sizeof(int));
      for (int i = 0; i < sizeV; i++)
        outFile.write((char*) & (*(svref[cpt]))(i), sizeof(double));
    }
    res = true;
  }
  else if (mode == "ascii")
  {
    outFile << nbVectors << endl;

    int sizeV = 0;
    for (int cpt = 0; cpt < nbVectors; cpt++)
    {
      sizeV = svref[cpt]->size();
      outFile << sizeV << endl;
      for (int i = 0; i < sizeV; i++)
      {
        /* WARNING this buffer is dangerous, the size is machine dependant*/
        char buffer[30];
        sprintf(buffer, "%1.17e ", (*(svref[cpt]))(i));
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


double* CompositeVector::getArray() const
{
  cout << "The SiconosVector is Composite, this function is not available" << endl;
  return NULL;
}


/*******************************************************************************
 *          GENERIC INTERNAL OPERATORS                                 *
 *******************************************************************************/

CompositeVector &CompositeVector::operator+=(const SiconosVector &v)
{
  //IN(" CompositeVector::operator+=(const SiconosVector &) \n");
  //cout<<" CompositeVector::operator+=(const SiconosVector &) \n";

  const int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" CompositeVector::operator+=   -- the vectors have not the same size");

  for (int i = 0; i < sizeV; i++)
    (*this)(i) += v(i);

  //OUT(" CompositeVector::operator+=(const SiconosVector &) \n");
  return *this;
}


CompositeVector &CompositeVector::operator-=(const SiconosVector &v)
{
  //IN(" CompositeVector::operator-=(const SiconosVector &) \n");
  //cout<<" CompositeVector::operator-=(const SiconosVector &) \n";

  const int sizeV = size();
  if (sizeV != v.size())
    SiconosVectorException::selfThrow(" CompositeVector::operator-=   -- the vectors have not the same size");

  for (int i = 0; i < sizeV; i++)
    (*this)(i) -= v(i);

  //OUT(" CompositeVector::operator-=(const SiconosVector &) \n");
  return *this;
}


CompositeVector& CompositeVector::operator = (const SiconosVector& v)
{
  //IN("CompositeVector::operator = \n");

  if (this != &v)
  {
    const int sizeV = size();
    if (sizeV != v.size())
      SiconosVectorException::selfThrow(" CompositeVector::operator = GENERIC  -- the vectors have not the same size");
    else
    {
      for (int i = 0; i < sizeV; i++)
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
    const int sizeV = size();
    if (sizeV == v.size())
    {
      if (sizeV == 0)
        res = true;
      else for (i = 0; i < sizeV; i++)
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
      const int sizeV = size();
      if (sizeV != v.size())
        SiconosVectorException::selfThrow(" CompositeVector::operator = GENERIC  -- the vectors have not the same size");
      else
      {
        for (int i = 0; i < sizeV; i++)
          (*this)(i) = v(i);
      }
    }
  }

  //OUT("CompositeVector::operator = \n");
  return *this;
}


/*******************************************************************************
 *          GENERIC EXTERNAL OPERATORS                                 *
/******************************************************************************/


/*******************************************************************************
 *          GENERIC EXTERNAL OPERATORS                                 *
/******************************************************************************/

CompositeVector CompositeVector::addition(const SiconosVector& v1) const
{
  //IN("CompositeVector::addition(const CompositeVector& v1)\n");

  const int sizeV = size();
  const int nbVec1 = svref.size();
  int i = 0;

  if (sizeV != v1.size())
    SiconosVectorException::selfThrow(" CompositeVector::addition(const CompositeVector& v1)  --- the vectors have not the same size");

  CompositeVector cv(*this);

  // the vectors have not the same structure
  cout << "WARNING : CompositeVector::addition(const CompositeVector& v1)  --- the vectors have not the same type (loops)" << endl;
  for (i = 0; i < sizeV; i++)
    cv(i) += v1(i);

  //OUT("CompositeVector::addition(const CompositeVector& v1)\n");
  return cv;
}


CompositeVector CompositeVector::subtraction(const SiconosVector& v1) const
{
  //IN("CompositeVector::subtraction(const CompositeVector& v1)\n");

  const int sizeV = size();
  const int nbVec1 = svref.size();
  int i = 0;

  if (sizeV != v1.size())
    SiconosVectorException::selfThrow(" CompositeVector::subtraction(const CompositeVector& v1)  --- the vectors have not the same size");

  CompositeVector cv(*this);

  // the vectors have not the same structure
  cout << "WARNING : CompositeVector::subtraction(const CompositeVector& v1)  --- the vectors have not the same type (loops)" << endl;
  for (i = 0; i < sizeV; i++)
    cv(i) -= v1(i);

  //OUT("CompositeVector::subtraction(const CompositeVector& v1)\n");
  return cv;
}

/*******************************************************************************
 *          SPECIFIC EXTERNAL OPERATORS                                 *
/******************************************************************************/

CompositeVector operator * (const CompositeVector& v, const double d)
{
  //IN(" CompositeVector operator * (const CompositeVector& v, const double d) \n");

  const int sizeV = v.svref.size();
  CompositeVector cv(v);
  for (int i = 0; i < sizeV; i++)
    *(v.svref[i]) *= d;

  //OUT(" CompositeVector operator * (const CompositeVector& v, const double d)\n");
  return cv;
}


CompositeVector operator * (const double d, const CompositeVector& v)
{
  //IN(" CompositeVector friend operator * \n");

  const int sizeV = v.svref.size();
  CompositeVector cv(v);
  for (int i = 0; i < sizeV; i++)
    *(v.svref[i]) *= d;

  //OUT(" CompositeVector friend operator * \n");
  return cv;
}


CompositeVector operator / (const CompositeVector& v, const double d)
{
  //IN(" CompositeVector friend operator * \n");

  if (d == 0.0)
    SiconosVectorException::selfThrow(" CompositeVector operator/   --- division by 0");

  const int sizeV = v.svref.size();
  CompositeVector cv(v);
  for (int i = 0; i < sizeV; i++)
    *(v.svref[i]) /= d;

  //OUT(" CompositeVector friend operator * \n");
  return cv;
}


CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2)
{
  //IN("CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2) \n");

  const int sizeV = v1.size();
  const int nbVec1 = v1.tabindex.size();
  const int nbVec2 = v2.tabindex.size();
  int i = 0;
  bool useBlas = true;

  if (sizeV != v2.size())
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
    for (i = 0; i < sizeV; i++)
      cv(i) += v2(i);
  }

  //OUT("CompositeVector operator + (const CompositeVector& v1, const CompositeVector& v2) \n");
  return cv;
}


CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2)
{
  //IN("CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2) \n");

  const int sizeV = v1.size();
  const int nbVec1 = v1.tabindex.size();
  const int nbVec2 = v2.tabindex.size();
  int i = 0;
  bool useBlas = true;

  if (sizeV != v2.size())
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
    for (i = 0; i < sizeV; i++)
      cv(i) -= v2(i);
  }

  //OUT("CompositeVector operator - (const CompositeVector& v1, const CompositeVector& v2) \n");
  return cv;
}

SimpleVector operator * (const SiconosMatrix &m, const CompositeVector &v)
{
  SimpleVector sv(m.size(0));

  LaVectorDouble lvd(v.size()), lvd2(sv.size());

  const int sizeV = v.size();

  for (int i = 0; i < sizeV; i++)
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

  const int sizeV = m.size(1);

  SimpleVector sv(sizeV);

  int i;
  LaVectorDouble lvd(v.size());
  LaVectorDouble lvd2(sizeV);

  for (i = 0; i < v.size(); i++)
    lvd(i) = v(i);
  Blas_Mat_Trans_Vec_Mult(m.getLaGenMatDouble(), lvd, lvd2, 1.0, 0.0);
  for (i = 0; i < sizeV; i++)
    sv(i) = lvd2(i);

  return sv;
}


