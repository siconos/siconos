
#include "TimeDiscretisation.h"

#include "check.h"

TimeDiscretisation::TimeDiscretisation()
{
  this->timeDiscretisationXML = NULL;
  this->k = 0;
  this->N = 0;
  this->h = 0.0;
  this->strategy = NULL;
  this->constant = NULL;
}

TimeDiscretisation::TimeDiscretisation(TimeDiscretisationXML *tdxml)
{
  this->timeDiscretisationXML = tdxml;
  this->strategy = NULL;
  this->k = 0;
}

TimeDiscretisation::~TimeDiscretisation()
{ }



void TimeDiscretisation::init(double t0, double T)
{
  IN("TimeDiscretisation::init\n");
  /* giving for example one of the following triplet :
   *    -  t0,T,h  --> Compute N
   *    -  t0,T,N --> Compute h
   *    -  t0,N,h --> Compute T
   */
  if ((t0 != -1) && (T != -1) && (this->h > 0))
  {
    this->N = floor((T - t0) / (this->h));
  }
  else if ((t0 != -1) && (T != -1) && (this->N > 0))
  {

  }
  else if ((t0 != -1) && (this->N != -1) && (this->h > 0))
  {

  }
  else
    RuntimeException::selfThrow("TimeDiscretisation::init(to, T) - insufficient data to init completely the time manager object");
  OUT("TimeDiscretisation::init\n");
}

void TimeDiscretisation::fillTimeDiscretisationWithTimeDiscretisationXML()
{
  IN("TimeDiscretisation::fillTimeDiscretisationWithTimeDiscretisationXML\n");
  if (this->timeDiscretisationXML != NULL)
  {
    /* \todo : check the coherency of the data */
    if (this->timeDiscretisationXML->hasH()) this->h = this->timeDiscretisationXML->getH();
    if (this->timeDiscretisationXML->hasN()) this->N = this->timeDiscretisationXML->getN();
    if (this->timeDiscretisationXML->hasTk())
    {
      this->tk = SimpleVector::SimpleVector();
      this->tk = *(this->timeDiscretisationXML->getTk());
    }
    if (this->timeDiscretisationXML->hasHMin()) this->hMin = this->timeDiscretisationXML->getHMin();
    if (this->timeDiscretisationXML->hasHMax()) this->hMax = this->timeDiscretisationXML->getHMax();
    this->constant = this->timeDiscretisationXML->isConstant();
  }
  else RuntimeException::selfThrow("TimeDiscretisation::fillTimeDiscretisationWithTimeDiscretisationXML - TimeDiscretisationXML object not exists");
  OUT("TimeDiscretisation::fillTimeDiscretisationWithTimeDiscretisationXML\n");
}

void TimeDiscretisation::display() const
{
  cout << "-----------------------------------------------------" << endl;
  cout << "____ data of the TimeDiscretisation " << endl;
  cout << "| h : " << this->h << endl;
  cout << "| N : " << this->N << endl;
  //cout<<"| tk size : "<<this->tk->size() <<" - "<<*(this->tk) <<endl;
  cout << "| tk " << endl;
  this->tk.display();
  cout << "| hMin : " << this->hMin << endl;
  cout << "| hMax : " << this->hMax << endl;
  cout << "| constant : " << this->constant << endl;
  cout << "-----------------------------------------------------" << endl << endl;
}

void TimeDiscretisation::saveTimeDiscretisationToXML()
{
  IN("TimeDiscretisation::saveTimeDiscretisationToXML\n");
  if (this->timeDiscretisationXML != NULL)
  {
    this->timeDiscretisationXML->setH(this->h);
    this->timeDiscretisationXML->setN(this->N);
    this->timeDiscretisationXML->setTk(&(this->tk));
    this->timeDiscretisationXML->setHMin(this->hMin);
    this->timeDiscretisationXML->setHMax(this->hMax);
  }
  else RuntimeException::selfThrow("TimeDiscretisation::saveTimeDiscretisationToXML - TimeDiscretisationXML object not exists");
  OUT("TimeDiscretisation::saveTimeDiscretisationToXML\n");
}

void TimeDiscretisation::createTimeDiscretisation(TimeDiscretisationXML * tdXML, double h, int N,
    SimpleVector *tk, double hMin, double hMax, bool constant, Strategy* str)//, OneStepIntegrator* osi)
//void TimeDiscretisation::createTimeDiscretisation(TimeDiscretisationXML * tdXML, Strategy* str)
{
  if (tdXML != NULL)
  {
    this->timeDiscretisationXML = tdXML;
    //  this->strategy = str;
    this->k = 0;

    this->fillTimeDiscretisationWithTimeDiscretisationXML();
    if (this->strategy != NULL) this->checkTimeDiscretisation();
    else RuntimeException::selfThrow("TimeDiscretisation::createTimeDiscretisation - no Strategy is defined for the TimeDiscretisation !");
  }
  else
  {
    if (str == NULL)
      RuntimeException::selfThrow("TimeDiscretisation::createTimeDiscretisation - The TimeDiscretisation belongs to no Strategy and no OneStepIntegrator !");
    else
    {
      /*
       * there's no required attributes to build the TimeDiscretisation
       */

      this->k = 0;
      this->h = h;
      this->N = N;

      this->tk = *tk;

      this->hMin = hMin;
      this->hMax = hMax;
      this->constant = constant;


      /*
       * to complete the TimeDiscretisation, we must know the node of the DOM tree
       * where we can build it
       */
      if (str != NULL)
      {
        this->strategy = str;
      }
      //      else
      //      {
      //        /*
      //         * \todo : the problem is for Mutli-Step Integrators
      //         */
      //      }
    }
  }
}

void TimeDiscretisation::checkTimeDiscretisation()
{
  IN("TimeDiscretisation::checkTimeDiscretisation\n");
  /* giving for example one of the following triplet :
   *    -  t0,T,h --> Compute N
   *    -  t0,T,N --> Compute h
   *    -  t0,N,h --> Compute T
   */

  double t0, T;

  t0 = this->strategy->getModel()->getT0();
  T  = this->strategy->getModel()->getFinalT();

  // T, N and h are optional attributes, whereas t0 is required
  bool hasT = false, hasT0 = true, hasN = false, hasH = false;

  hasT = this->strategy->getModel()->getSiconosModelXML()->hasT();
  hasN = this->timeDiscretisationXML->hasN();
  hasH = this->timeDiscretisationXML->hasH();

  if (!hasT0 || !hasT || !hasN || !hasH)
  {
    if (hasT0 && hasT && hasH)
    {
      /*
       *  N must be computed
       */
      this->N = floor(T - t0) / (this->h);
      //      cout<<"N == "<<floor(T - t0)/(this->h)<<"                press ENTER"<<endl;getchar();
    }
    else if (hasT0 && hasT && hasN)
    {
      /*
       *  h must be computed
       */
      this->h = floor(T - t0) / ((double)N);
      //      cout<<"h == "<<floor(T-t0) / ((double)N)<<"                press ENTER"<<endl;getchar();
    }
    else if (hasT0 && hasN && hasH)
    {
      /*
       *  T must be computed
       */
      this->strategy->getModel()->setFinalT(t0 + ((double)this->N) * (this->h));
      //      cout<<"T == "<<t0 + ((double)this->N) * (this->h)<<"                press ENTER"<<endl;getchar();
    }
    else
      RuntimeException::selfThrow("TimeDiscretisation::checkTimeDiscretisation - insufficient data to init completely the time manager object, you need : (t0,T,h), (t0,T,N) or (t0,h,N)");
  }
  // \todo : void TimeDiscretisation::checkTimeDiscretisation() - checking the coherency of the data given
  OUT("TimeDiscretisation::checkTimeDiscretisation\n");
}
