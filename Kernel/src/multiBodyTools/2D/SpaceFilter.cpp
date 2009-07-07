#include "SpaceFilter.hpp"


/* hash is done with encapsulation */

class Hashed : public boost::enable_shared_from_this<Hashed>
{
public:
  SP::CircularDS body;
  int i;
  int j;
  int k;
  Hashed(SP::CircularDS body, int i, int j, int k = 0) :
    body(body), i(i), j(j), k(k) {};

  ~Hashed() {};

};

class HashedCircle : public Hashed
{
public:
  SP::Circle body;
  HashedCircle(SP::Circle c, int i, int j) : Hashed(c, i, j) {};

};

class HashedDisk : public Hashed
{
public:
  SP::Disk body;
  HashedDisk(SP::Disk d, int i, int j) : Hashed(d, i, j) {};

};

/* needed by boost hash */
bool operator ==(SP::Hashed const& a, SP::Hashed const& b)
{
  return (a->i == b->i &&
          a->j == b->j &&
          a->k == b->k);
};



std::size_t hash_value(SP::Hashed const& h)
{
  std::size_t seed = 0;
  boost::hash_combine(seed, h->i);
  boost::hash_combine(seed, h->j);
  boost::hash_combine(seed, h->k);
  return seed;
};


TYPEDEF_SPTR(HashedCircle);
TYPEDEF_SPTR(HashedDisk);

/* the hashing is done with a visitor */
class SpaceFilter::_BodyHash : public SiconosVisitor
{
public:
  SP::SpaceFilter parent;
  _BodyHash(SP::SpaceFilter p) : parent(p) {};

  void visit(SP::Disk pds)
  {
    int i, j, imin, imax, jmin, jmax;

    unsigned int _bboxfactor = parent->bboxfactor();
    unsigned int _cellsize = parent->cellsize();

    imin = (int) floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int) floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int) floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int) floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i)
    {
      for (j = jmin; j <= jmax; ++j)
      {
        parent->insert(pds, i, j, 0);
      }
    }
  };

  void visit(SP::Circle pds)
  {
    int i, j, imin, imax, jmin, jmax;

    unsigned int _bboxfactor = parent->bboxfactor();
    unsigned int _cellsize = parent->cellsize();

    imin = (int) floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int) floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int) floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int) floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i)
    {
      for (j = jmin; j <= jmax; ++j)
      {
        parent->insert(pds, i, j, 0);
      }
    }
  }

  /* ... visit other objects*/

};




/* proximity detection for circular object */
void SpaceFilter::_CircularCircularFilter(SP::CircularDS ds1,
    SP::CircularDS ds2)
{

  SP::CircularR rel;
  SP::DynamicalSystemsGraph DSG0 = _nsds->getTopologyPtr()->getDSGPtr(0);

  assert(ds1 != ds2);
  assert(DSG0->bundle(DSG0->descriptor(ds1)) == ds1);
  assert(DSG0->bundle(DSG0->descriptor(ds2)) == ds2);

  double r1 = ds1->getRadius();
  double r2 = ds2->getRadius();
  double tol = r1 + r2;

  double x1 = ds1->getQ(0);
  double y1 = ds1->getQ(1);
  double x2 = ds2->getQ(0);
  double y2 = ds2->getQ(1);
  double rmax = fmax(r1, r2);
  double rmin = fmin(r1, r2);

  double d = hypot(x1 - x2, y1 - y2);

  if (d < rmax)
  {
    // one inside the other
    if (rmax - (d + rmin) < tol)
      rel.reset(new CircleCircleR(r1, r2));
  }
  else
  {
    if (d - (r1 + r2) < tol)
      rel.reset(new DiskDiskR(r1, r2));
  }


  if (rel)
  {
    bool found = false;
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (boost::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
         oei != oeiend; ++oei)
    {
      if (DSG0->bundle(DSG0->target(*oei)) == ds2)
      {
        found = true;
        break;
      }
    }

    if (!found)
    {
      SP::Interaction inter(new Interaction(_interID++,
                                            2,
                                            _nslaw,
                                            rel));
      inter->insert(ds1);
      inter->insert(ds2);
      _nsds->getTopologyPtr()->addInteraction(inter);
    }
  }
  else
  {
    // is interaction in graph ?
    bool found = false;
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (boost::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
         oei != oeiend; ++oei)
    {
      if (DSG0->bundle(DSG0->target(*oei)) == ds2)
      {
        found = true;
        break;
      }
    }

    if (found)
    {
      _nsds->getTopologyPtr()->
      removeInteraction(DSG0->bundle(*oei)->getInteractionPtr());
    }

  }
};

/* disk plan relation comparison */
struct SpaceFilter::_IsSameDiskPlanR : public SiconosVisitor
{
  double A, B, C, r, xCenter, yCenter, width;
  bool flag;
  SP::SpaceFilter parent;
  _IsSameDiskPlanR(SP::SpaceFilter p, double A, double B, double C, double r,
                   double xCenter, double yCenter, double width) :
    parent(p), A(A), B(B), C(C), r(r), xCenter(xCenter), yCenter(yCenter), width(width), flag(false) {};


  void visit(SP::DiskDiskR)
  {
    flag = false;
  };

  void visit(SP::CircleCircleR)
  {
    flag = false;
  };

  void visit(SP::DiskPlanR rel)
  {
    flag = rel->equal(A, B, C, r, xCenter, yCenter, width);
  };


};


/* proximity detection between circular object and plans */
void SpaceFilter::_PlanCircularFilter(double A, double B, double C,
                                      double xCenter, double yCenter, double width, SP::CircularDS ds)
{
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  SP::DynamicalSystemsGraph DSG0 = _nsds->getTopologyPtr()->getDSGPtr(0);

  boost::shared_ptr<_IsSameDiskPlanR> isSameDiskPlanR(new _IsSameDiskPlanR(shared_from_this(), A, B, C, r, xCenter, yCenter, width));


  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  SP::DiskPlanR relp(new DiskPlanR(r, A, B, C, xCenter, yCenter, width));
  if (relp->distance(ds->getQ(0),
                     ds->getQ(1),
                     ds->getRadius()) < tol)

  {
    // is interaction in graph ?
    bool found = false;
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (boost::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)->getInteractionPtr()
      ->getRelationPtr()->accept(isSameDiskPlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameDiskPlanR->flag)
      {
        found = true;
        break;
      }
    }
    if (!found)
      // no
    {
      SP::Interaction inter(new Interaction(_interID++,
                                            2,
                                            _nslaw,
                                            relp));
      inter->insert(ds);

      _nsds->getTopologyPtr()->addInteraction(inter);

    }
  }
  else
  {
    // is interaction in graph ?
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (boost::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)->getInteractionPtr()
      ->getRelationPtr()->accept(isSameDiskPlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameDiskPlanR->flag)
      {
        _nsds->getTopologyPtr()->
        removeInteraction(DSG0->bundle(*oei)->getInteractionPtr());
        break;
      }
    }
  }
};


/* insertion */
void SpaceFilter::insert(SP::Disk ds, int i, int j, int k)
{

  SP::Hashed hashed(new HashedDisk(ds, i, j));
  _hash_table.insert(hashed);
};

void SpaceFilter::insert(SP::Circle ds, int i, int j, int k)
{

  SP::Hashed hashed(new HashedCircle(ds, i, j));
  _hash_table.insert(hashed);
};

/* insert other objects */



/* dynamical systems proximity detection */
typedef std::pair<int, int> interPair;

bool operator ==(interPair const& a, interPair const& b)
{
  return ((a.first == b.first) && (a.second == b.second));
};

struct SpaceFilter::_FindInteractions : public SiconosVisitor
{

  typedef std::tr1::unordered_multiset<interPair, boost::hash<interPair> > interPairs;

  SP::SpaceFilter parent;
  _FindInteractions(SP::SpaceFilter p) : parent(p) {};

  void visit_circular(SP::CircularDS  ds1)
  {
    assert(parent->_plans->size(0) > 0);
    for (unsigned int i = 0; i < parent->_plans->size(0); ++i)
    {

      // to avoid cast see dual dispatch, visitor pattern
      parent->_PlanCircularFilter((*parent->_plans)(i, 0), (*parent->_plans)(i, 1), (*parent->_plans)(i, 2),
                                  (*parent->_plans)(i, 3), (*parent->_plans)(i, 4), (*parent->_plans)(i, 5),
                                  ds1);
    }

    SP::SiconosVector Q1 = ds1->getQPtr();

    double x1 = Q1->getValue(0);
    double y1 = Q1->getValue(1);
    SP::Hashed hds1(new Hashed(ds1, (int) floor(x1 / parent->_cellsize),
                               (int) floor(y1 / parent->_cellsize)));
    std::pair<space_hash::iterator, space_hash::iterator>
    neighbours = parent->_hash_table.equal_range(hds1);

    unsigned int j;
    interPairs declaredInteractions;

    for (j = 0; neighbours.first != neighbours.second; ++neighbours.first, ++j)
    {
      SP::CircularDS ds2 = (*neighbours.first)->body;
      int ids1 = ds1->getNumber();
      int ids2 = ds2->getNumber();
      int imax = (std::max)(ids1, ids2);
      int imin = (std::min)(ids1, ids2);
      if (ids1 != ids2)
      {
        interPair interpair;

        interpair = std::pair<int, int>(imin, imax);
        if (declaredInteractions.find(interpair)
            == declaredInteractions.end())
        {

          declaredInteractions.insert(interpair);

          parent->_CircularCircularFilter(ds1, ds2);
        }

      }
    }

  };

  void visit(SP::Circle circle)
  {
    visit_circular(circle);
  };


  void visit(SP::Disk disk)
  {
    visit_circular(disk);
  };
};


/* general proximity detection */
void SpaceFilter::buildInteractions()
{

  DSIterator it1;
  unsigned int i, j;


  SP::DynamicalSystemsGraph DSG0 = _nsds->getTopologyPtr()->getDSGPtr(0);

  boost::shared_ptr<_BodyHash> hasher(new _BodyHash(shared_from_this()));
  boost::shared_ptr<_FindInteractions> findInteractions(new _FindInteractions(shared_from_this()));


  _hash_table.clear();

  // 1: rehash DS
  DynamicalSystemsGraph::VIterator vi, viend;
  for (boost::tie(vi, viend) = DSG0->vertices();
       vi != viend; ++vi)
  {
    // to avoid cast see dual dispatch, visitor pattern
    DSG0->bundle(*vi)->accept(hasher);
  }

  // 2: prox detection
  for (boost::tie(vi, viend) = DSG0->vertices();
       vi != viend; ++vi)
  {

    DSG0->bundle(*vi)->accept(findInteractions);

  }

}
