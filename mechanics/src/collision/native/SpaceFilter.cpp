#include <boost/numeric/bindings/ublas/matrix.hpp>

#include "SpaceFilter.hpp"
#include "SpaceFilter_impl.hpp"

#include "Circle.hpp"
#include "Disk.hpp"
#include "DiskDiskR.hpp"
#include "CircleCircleR.hpp"
#include "DiskPlanR.hpp"
#include "DiskMovingPlanR.hpp"
#include "SphereLDS.hpp"
#include "SphereLDSSphereLDSR.hpp"
#include "SphereNEDSSphereNEDSR.hpp"
#include "SphereLDSPlanR.hpp"
#include "SphereNEDS.hpp"
#include "SphereNEDSPlanR.hpp"
#include "ExternalBody.hpp"

#include <Model.hpp>
#include <Simulation.hpp>
#include <NonSmoothDynamicalSystem.hpp>
#include <SimulationTypeDef.hpp>
#include <NonSmoothLaw.hpp>


#include <cmath>
//#define DEBUG_MESSAGES 1
#include "debug.h"





/* hash is done with encapsulation */


/* needed by boost hash */
bool operator ==(SP::Hashed const& a, SP::Hashed const& b);
bool operator ==(SP::Hashed const& a, SP::Hashed const& b)
{
  return (a->i == b->i &&
          a->j == b->j &&
          a->k == b->k);
}


std::size_t hash_value(SP::Hashed const& h);
std::size_t hash_value(SP::Hashed const& h)
{
  std::size_t seed = 0;
  boost::hash_combine(seed, h->i);
  boost::hash_combine(seed, h->j);
  boost::hash_combine(seed, h->k);
  return seed;
}

SpaceFilter::SpaceFilter(unsigned int bboxfactor,
              unsigned int cellsize,
              SP::Model model,
              SP::SiconosMatrix plans,
              SP::FMatrix moving_plans) :
    _bboxfactor(bboxfactor),
    _cellsize(cellsize),
    _model(model),
    _nslaws(new NSLawMatrix()),
    _plans(plans),
    _moving_plans(moving_plans),
    _osnsinit(false),
    _hash_table(new space_hash()),
    diskdisk_relations(new DiskDiskRDeclaredPool()),
    diskplan_relations(new DiskPlanRDeclaredPool()),
  circlecircle_relations(new CircleCircleRDeclaredPool())
{};

SpaceFilter::SpaceFilter(unsigned int bboxfactor,
              unsigned int cellsize,
              SP::Model model,
              SP::SiconosMatrix plans) :
    _bboxfactor(bboxfactor),
    _cellsize(cellsize),
    _model(model),
    _nslaws(new NSLawMatrix()),
    _plans(plans),
    _osnsinit(false),
    _hash_table(new space_hash()),
    diskdisk_relations(new DiskDiskRDeclaredPool()),
    diskplan_relations(new DiskPlanRDeclaredPool()),
    circlecircle_relations(new CircleCircleRDeclaredPool())
{};

SpaceFilter::SpaceFilter(SP::Model model) :
  _model(model),
  _nslaws(new NSLawMatrix()),
  _osnsinit(false),
  _hash_table(new space_hash()),
  diskdisk_relations(new DiskDiskRDeclaredPool()),
  diskplan_relations(new DiskPlanRDeclaredPool()),
  circlecircle_relations(new CircleCircleRDeclaredPool())
{};

SpaceFilter::SpaceFilter() :
  _osnsinit(false),
  _hash_table(new space_hash()),
  diskdisk_relations(new DiskDiskRDeclaredPool()),
  diskplan_relations(new DiskPlanRDeclaredPool()),
  circlecircle_relations(new CircleCircleRDeclaredPool())
{};


/* the hashing is done with a visitor */
struct SpaceFilter::_BodyHash : public SiconosVisitor
{
public:
  SpaceFilter& parent;
  _BodyHash(SpaceFilter& p) : parent(p) {};

  using SiconosVisitor::visit;

  void visit(SP::Disk pds)
  {
    int i, j, imin, imax, jmin, jmax;
    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int) floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int) floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int) floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int) floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i)
    {
      for (j = jmin; j <= jmax; ++j)
      {
        parent.insert(pds, i, j, 0);
      }
    }
  };

  void visit(SP::Circle pds)
  {
    int i, j, imin, imax, jmin, jmax;

    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int) floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int) floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int) floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int) floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i)
    {
      for (j = jmin; j <= jmax; ++j)
      {
        parent.insert(pds, i, j, 0);
      }
    }
  }

  void visit(SP::SphereLDS pds)
  {
    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;

    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int) floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int) floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int) floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int) floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    kmin = (int) floor((pds->getQ(2) - _bboxfactor * pds->getRadius()) / _cellsize);
    kmax = (int) floor((pds->getQ(2) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i)
    {
      for (j = jmin; j <= jmax; ++j)
      {
        for (k = kmin; k <= kmax; ++k)
        {
          parent.insert(pds, i, j, k);
        }
      }
    }

  }

  void visit(SP::SphereNEDS pds)
  {
    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;

    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int) floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int) floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int) floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int) floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    kmin = (int) floor((pds->getQ(2) - _bboxfactor * pds->getRadius()) / _cellsize);
    kmax = (int) floor((pds->getQ(2) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i)
    {
      for (j = jmin; j <= jmax; ++j)
      {
        for (k = kmin; k <= kmax; ++k)
        {
          parent.insert(pds, i, j, k);
        }
      }
    }

  }

  void visit(SP::ExternalBody d)
  {
    d->selfHash(parent);
  }

  /* ... visit other objects */

};



/* proximity detection for circular object */
struct SpaceFilter::_CircularFilter : public SiconosVisitor
{
  using SiconosVisitor::visit;

  SP::SpaceFilter parent;
  SP::CircularDS ds1;

  _CircularFilter(SP::SpaceFilter parent, SP::CircularDS ds1) : parent(parent), ds1(ds1) {};

  void visit_circular(SP::CircularDS ds2)
  {

    SP::CircularR rel;
    SP::DynamicalSystemsGraph DSG0 = parent->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

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
      // one inside the other : CircleCircle relation
      if (rmax - (d + rmin) < tol)
      {
        CircleCircleRDeclaredPool::iterator rcandid =
          parent->circlecircle_relations->find(CircleCircleRDeclared(r1, r2));
        if (rcandid == parent->circlecircle_relations->end())
        {
          // a new relation
          rel.reset(new CircleCircleR(r1, r2));
          (*(parent->circlecircle_relations))[CircleCircleRDeclared(r1, r2)] = rel;
        }
        else
        {
          // get relation from pool
          rel = (*rcandid).second;
        }
      }
    }
    else
    {
      // a DiskDisk relation
      if (d - (r1 + r2) < tol)
      {
        DiskDiskRDeclaredPool::iterator rcandid =
          parent->diskdisk_relations->find(DiskDiskRDeclared(r1, r2));
        if (rcandid == parent->diskdisk_relations->end())
        {
          // a new relation
          rel.reset(new DiskDiskR(r1, r2));

          // FIX : this does not work!!
          // parent->diskdisk_relations[DiskDiskRDeclared(r1,r2)] = rel;
        }
        else
        {
          // get relation from pool
          rel = (*rcandid).second;
        }
      }
    }

    if (rel)
    {
      bool found = false;
      DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
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
        SP::NonSmoothLaw nslaw = (*parent->_nslaws)(DSG0->groupId[DSG0->descriptor(ds1)],
                                                    DSG0->groupId[DSG0->descriptor(ds2)]);

        SP::Interaction inter(new Interaction(nslaw, rel));
        parent->model()->simulation()->link(inter, ds1, ds2);
      }
    }
    else
    {
      // is interaction in graph ?
      bool found = false;
      DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
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


        DEBUG_PRINTF("remove interaction : %d\n", DSG0->bundle(*oei)->number());
        parent->model()->nonSmoothDynamicalSystem()->topology()->
        removeInteraction(DSG0->bundle(*oei));
      }
    }
  }

  void visit(SP::Disk disk)
  {
    visit_circular(disk);
  }

  void visit(SP::Circle circle)
  {
    visit_circular(circle);
  }

  // do nothing (everything must be done in ExternalBody.findInteractions)
  void visit(SP::ExternalBody)
  {}


};

/* proximity detection for sphere objects */
struct SpaceFilter::_SphereLDSFilter : public SiconosVisitor
{
  using SiconosVisitor::visit;

  SP::SpaceFilter parent;
  SP::SphereLDS ds1;

  _SphereLDSFilter(SP::SpaceFilter p, SP::SphereLDS s) : parent(p), ds1(s) {};

  void visit(SP::SphereLDS ds2)
  {
    SP::SphereLDSSphereLDSR rel;
    SP::DynamicalSystemsGraph DSG0 = parent->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

    assert(ds1 != ds2);
    assert(DSG0->bundle(DSG0->descriptor(ds1)) == ds1);
    assert(DSG0->bundle(DSG0->descriptor(ds2)) == ds2);

    double r1 = ds1->getRadius();
    double r2 = ds2->getRadius();
    double tol = r1 + r2;

    double x1 = ds1->getQ(0);
    double y1 = ds1->getQ(1);
    double z1 = ds1->getQ(2);
    double x2 = ds2->getQ(0);
    double y2 = ds2->getQ(1);
    double z2 = ds2->getQ(2);

    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;

    double d = sqrt(dx * dx + dy * dy + dz * dz);

    if (d < 2 * tol)
    {
      rel.reset(new SphereLDSSphereLDSR(r1, r2));

      bool found = false;
      DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
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
        SP::NonSmoothLaw nslaw = (*parent->_nslaws)(DSG0->groupId[DSG0->descriptor(ds1)],
                                                    DSG0->groupId[DSG0->descriptor(ds2)]);

        SP::Interaction inter(new Interaction(nslaw, rel));
        parent->model()->simulation()->link(inter, ds1, ds2);
      }
    }
    else
    {
      // is interaction in graph ?
      bool found = false;
      DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
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
        parent->model()->nonSmoothDynamicalSystem()->topology()->
        removeInteraction(DSG0->bundle(*oei));
      }

    }
  }
};


struct SpaceFilter::_SphereNEDSFilter : public SiconosVisitor
{
  using SiconosVisitor::visit;

  SP::SpaceFilter parent;
  SP::SphereNEDS ds1;

  _SphereNEDSFilter(SP::SpaceFilter p, SP::SphereNEDS s) : parent(p), ds1(s) {};

  void visit(SP::SphereNEDS ds2)
  {
    SP::SphereNEDSSphereNEDSR rel;
    SP::DynamicalSystemsGraph DSG0 = parent->model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

    assert(ds1 != ds2);
    assert(DSG0->bundle(DSG0->descriptor(ds1)) == ds1);
    assert(DSG0->bundle(DSG0->descriptor(ds2)) == ds2);

    double r1 = ds1->getRadius();
    double r2 = ds2->getRadius();
    double tol = r1 + r2;

    double x1 = ds1->getQ(0);
    double y1 = ds1->getQ(1);
    double z1 = ds1->getQ(2);
    double x2 = ds2->getQ(0);
    double y2 = ds2->getQ(1);
    double z2 = ds2->getQ(2);

    double dx = x1 - x2;
    double dy = y1 - y2;
    double dz = z1 - z2;

    double d = sqrt(dx * dx + dy * dy + dz * dz);

    if (d < 2 * tol)
    {
      rel.reset(new SphereNEDSSphereNEDSR(r1, r2));

      bool found = false;
      DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
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
        SP::NonSmoothLaw nslaw = (*parent->_nslaws)(DSG0->groupId[DSG0->descriptor(ds1)],
                                                    DSG0->groupId[DSG0->descriptor(ds2)]);

        SP::Interaction inter(new Interaction(nslaw, rel));

        parent->model()->simulation()->link(inter, ds1, ds2);
      }
    }
    else
    {
      // is interaction in graph ?
      bool found = false;
      DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1));
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
        parent->model()->nonSmoothDynamicalSystem()->topology()->
        removeInteraction(DSG0->bundle(*oei));
      }

    }
  }
};



/* disk plan relation comparison */
struct SpaceFilter::_IsSameDiskPlanR : public SiconosVisitor
{

  using SiconosVisitor::visit;

  SP::SpaceFilter parent;
  double A, B, C, r, xCenter, yCenter, width;
  bool flag;
  _IsSameDiskPlanR(SP::SpaceFilter p, double A, double B, double C, double r,
                   double xCenter, double yCenter, double width) :
    parent(p), A(A), B(B), C(C), r(r), xCenter(xCenter), yCenter(yCenter), width(width), flag(false) {};


  void visit(const DiskDiskR&)
  {
    flag = false;
  };

  void visit(const CircleCircleR&)
  {
    flag = false;
  };

  void visit(const DiskMovingPlanR&)
  {
    flag = false;
  };

  void visit(const LagrangianScleronomousR&)
  {
    flag = false;
  };

  void visit(const DiskPlanR& rel)
  {
    flag = rel.equal(A, B, C, r, xCenter, yCenter, width);
  };

};

struct SpaceFilter::_IsSameDiskMovingPlanR : public SiconosVisitor
{

  using SiconosVisitor::visit;

  SP::SpaceFilter parent;
  FTime AF, BF, CF;
  double r;
  bool flag;
  _IsSameDiskMovingPlanR(SP::SpaceFilter p, FTime AF, FTime BF, FTime CF, double r) :
    parent(p), AF(AF), BF(BF), CF(CF), r(r), flag(false) {};


  void visit(const DiskDiskR&)
  {
    flag = false;
  };

  void visit(const CircleCircleR&)
  {
    flag = false;
  };

  void visit(const DiskPlanR&)
  {
    flag = false;
  };

  void visit(const LagrangianScleronomousR&)
  {
    flag = false;
  }

  void visit(const DiskMovingPlanR& rel)
  {
    flag = rel.equal(AF, BF, CF, r);
  };

};

/* sphere plan relation comparison */
struct SpaceFilter::_IsSameSpherePlanR : public SiconosVisitor
{

  using SiconosVisitor::visit;

  SP::SpaceFilter parent;
  double A, B, C, D, r;
  bool flag;
  _IsSameSpherePlanR(SP::SpaceFilter p, double A, double B, double C, double D, double r):
    parent(p), A(A), B(B), C(C), D(D), r(r), flag(false) {};

  void visit(const SphereLDSSphereLDSR&)
  {
    flag = false;
  };

  void visit(const SphereNEDSSphereNEDSR&)
  {
    flag = false;
  };

  void visit(const SphereLDSPlanR& rel)
  {
    flag = rel.equal(A, B, C, D, r);
  };

  void visit(const SphereNEDSPlanR& rel)
  {
    flag = rel.equal(A, B, C, D, r);
  };
};


/* proximity detection between circular object and plans */
void SpaceFilter::_PlanCircularFilter(double A, double B, double C,
                                      double xCenter, double yCenter,
                                      double width, SP::CircularDS ds)
{
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  SP::DynamicalSystemsGraph DSG0 = model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

  _IsSameDiskPlanR
  isSameDiskPlanR = _IsSameDiskPlanR(shared_from_this(), A, B, C, r,
                                     xCenter, yCenter, width);


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
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(isSameDiskPlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameDiskPlanR.flag)
      {
        found = true;
        break;
      }
    }

    if (!found)
      // no
    {
      SP::NonSmoothLaw nslaw = (*_nslaws)(DSG0->groupId[DSG0->descriptor(ds)],
                                          DSG0->groupId[DSG0->descriptor(ds)]);

      SP::Interaction inter(new Interaction(nslaw, relp));
      DEBUG_PRINTF("insert interaction : %d\n", inter->number());
      model()->simulation()->link(inter, ds);
    }
  }
  else
  {
    // is interaction in graph ?
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(isSameDiskPlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameDiskPlanR.flag)
      {

        DEBUG_PRINTF("remove interaction : %d\n", DSG0->bundle(*oei)->number());

        model()->nonSmoothDynamicalSystem()->topology()->
        removeInteraction(DSG0->bundle(*oei));
        break;
      }
    }
  }
}



/* proximity detection between circular object and plans */
void SpaceFilter::_MovingPlanCircularFilter(unsigned int i, SP::CircularDS ds, double time)
{
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  SP::DynamicalSystemsGraph DSG0 = model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

  _IsSameDiskMovingPlanR
  isSameDiskMovingPlanR = _IsSameDiskMovingPlanR(shared_from_this(),
                          (*_moving_plans)(i, 0),
                          (*_moving_plans)(i, 1),
                          (*_moving_plans)(i, 2),
                          r);

  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  SP::DiskMovingPlanR relp(new DiskMovingPlanR((*_moving_plans)(i, 0), (*_moving_plans)(i, 1), (*_moving_plans)(i, 2),
                           (*_moving_plans)(i, 3), (*_moving_plans)(i, 4), (*_moving_plans)(i, 5), r));

  relp->init(time);

  if (relp->distance(ds->getQ(0),
                     ds->getQ(1),
                     ds->getRadius()) < tol)

  {
    // is interaction in graph ?
    bool found = false;
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(isSameDiskMovingPlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameDiskMovingPlanR.flag)
      {
        found = true;
        break;
      }
    }
    if (!found)
      // no
    {
      SP::NonSmoothLaw nslaw = (*_nslaws)(DSG0->groupId[DSG0->descriptor(ds)],
                                          DSG0->groupId[DSG0->descriptor(ds)]);

      SP::Interaction inter(new Interaction(nslaw, relp));
      model()->simulation()->link(inter, ds);
    }
  }
  else
  {
    // is interaction in graph ?
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(isSameDiskMovingPlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameDiskMovingPlanR.flag)
      {
        model()->nonSmoothDynamicalSystem()->topology()->
        removeInteraction(DSG0->bundle(*oei));
        break;
      }
    }
  }
}

/* proximity detection between circular object and plans */
void SpaceFilter::_PlanSphereLDSFilter(double A, double B, double C, double D, SP::SphereLDS ds)
{
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  SP::DynamicalSystemsGraph DSG0 = model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

  _IsSameSpherePlanR
  IsSameSpherePlanR =
    _IsSameSpherePlanR(shared_from_this(), A, B, C, D, r);


  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  SP::SphereLDSPlanR relp(new SphereLDSPlanR(r, A, B, C, D));
  if (relp->distance(ds->getQ(0),
                     ds->getQ(1),
                     ds->getQ(2),
                     ds->getRadius()) < tol)

  {
    // is interaction in graph ?
    bool found = false;
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(IsSameSpherePlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && IsSameSpherePlanR.flag)
      {
        found = true;
        break;
      }
    }
    if (!found)
      // no
    {
      SP::NonSmoothLaw nslaw = (*_nslaws)(DSG0->groupId[DSG0->descriptor(ds)],
                                                  DSG0->groupId[DSG0->descriptor(ds)]);

      SP::Interaction inter(new Interaction(nslaw, relp));
      model()->simulation()->link(inter, ds);
    }
  }
  else
  {
    // is interaction in graph ?
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(IsSameSpherePlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && IsSameSpherePlanR.flag)
      {
        model()->nonSmoothDynamicalSystem()->topology()->
        removeInteraction(DSG0->bundle(*oei));
        break;
      }
    }
  }
}


// note : all PlanObject should be merged
void SpaceFilter::_PlanSphereNEDSFilter(double A, double B, double C, double D, SP::SphereNEDS ds)
{
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  SP::DynamicalSystemsGraph DSG0 = model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

  _IsSameSpherePlanR
  isSameSpherePlanR =
    _IsSameSpherePlanR(shared_from_this(), A, B, C, D, r);


  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  SP::SphereNEDSPlanR relp(new SphereNEDSPlanR(r, A, B, C, D));
  if (relp->distance(ds->getQ(0),
                     ds->getQ(1),
                     ds->getQ(2),
                     ds->getRadius()) < tol)

  {
    // is interaction in graph ?
    bool found = false;
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(isSameSpherePlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameSpherePlanR.flag)
      {
        found = true;
        break;
      }
    }
    if (!found)
      // no
    {
      SP::NonSmoothLaw nslaw = (*_nslaws)(DSG0->groupId[DSG0->descriptor(ds)],
                                          DSG0->groupId[DSG0->descriptor(ds)]);

      SP::Interaction inter(new Interaction(nslaw, relp));
      model()->simulation()->link(inter, ds);
    }
  }
  else
  {
    // is interaction in graph ?
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std11::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds));
         oei != oeiend; ++oei)
    {
      DSG0->bundle(*oei)
      ->relation()->accept(isSameSpherePlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds
          && isSameSpherePlanR.flag)
      {
        model()->nonSmoothDynamicalSystem()->topology()->
        removeInteraction(DSG0->bundle(*oei));
        break;
      }
    }
  }
}


/* insertion */
void SpaceFilter::insert(SP::Disk ds, int i, int j, int k)
{

  SP::Hashed hashed(new Hashed(std11::static_pointer_cast<DynamicalSystem>(ds), i, j));
  _hash_table->insert(hashed);
}

void SpaceFilter::insert(SP::Circle ds, int i, int j, int k)
{

  SP::Hashed hashed(new Hashed(std11::static_pointer_cast<DynamicalSystem>(ds), i, j));
  _hash_table->insert(hashed);
}

void SpaceFilter::insert(SP::SphereLDS ds, int i, int j, int k)
{

  SP::Hashed hashed(new Hashed(ds, i, j, k));
  _hash_table->insert(hashed);
}

void SpaceFilter::insert(SP::SphereNEDS ds, int i, int j, int k)
{

  SP::Hashed hashed(new Hashed(ds, i, j, k));
  _hash_table->insert(hashed);
}

void SpaceFilter::insert(SP::Hashed hashed)
{
  _hash_table->insert(hashed);
}

/* insert other objects */



/* dynamical systems proximity detection */
typedef std::pair<int, int> interPair;
bool operator ==(interPair const& a, interPair const& b);
bool operator ==(interPair const& a, interPair const& b)
{
  return ((a.first == b.first) && (a.second == b.second));
}

bool operator ==(std::pair<double, double> const& a,
                 std::pair<double, double> const& b);
bool operator ==(std::pair<double, double> const& a,
                 std::pair<double, double> const& b)
{
  return ((a.first == b.first) && (a.second == b.second));
}

bool operator ==(DiskPlanRDeclared const& a, DiskPlanRDeclared const& b);
bool operator ==(DiskPlanRDeclared const& a, DiskPlanRDeclared const& b)
{
  return ((a[0] == b[0] &&
           a[1] == b[1] &&
           a[2] == b[2] &&
           a[3] == b[3] &&
           a[4] == b[4] &&
           a[5] == b[5]));
}

struct SpaceFilter::_FindInteractions : public SiconosVisitor
{

  using SiconosVisitor::visit;

  typedef boost::unordered_multiset < interPair,
          boost::hash<interPair> > interPairs;

  SP::SpaceFilter parent;
  double time;
  _FindInteractions(SP::SpaceFilter p, double time) : parent(p), time(time) {};

  void visit_circular(SP::CircularDS  ds1)
  {
    assert(parent->_plans->size(0) > 0);

    // interactions with plans

    if (parent->_plans)
    {
      for (unsigned int i = 0; i < parent->_plans->size(0); ++i)
      {
        parent->_PlanCircularFilter((*parent->_plans)(i, 0),
                                    (*parent->_plans)(i, 1),
                                    (*parent->_plans)(i, 2),
                                    (*parent->_plans)(i, 3),
                                    (*parent->_plans)(i, 4),
                                    (*parent->_plans)(i, 5),
                                    ds1);
      }
    }

    if (parent->_moving_plans)
    {
      for (unsigned int i = 0; i < parent->_moving_plans->size1(); ++i)
      {
        parent->_MovingPlanCircularFilter(i, ds1, time);
      }
    }

    SP::SiconosVector Q1 = ds1->q();

    double x1 = Q1->getValue(0);
    double y1 = Q1->getValue(1);
    SP::Hashed hds1(new Hashed(std11::static_pointer_cast<DynamicalSystem>(ds1),
                             (int) floor(x1 / parent->_cellsize),
                               (int) floor(y1 / parent->_cellsize)));

    // find all other systems that are in the same cells
    std::pair<space_hash::iterator, space_hash::iterator>
    neighbours = parent->_hash_table->equal_range(hds1);

    unsigned int j;
    interPairs declaredInteractions;
    std11::shared_ptr<_CircularFilter>
    circularFilter(new _CircularFilter(parent, ds1));

    for (j = 0; neighbours.first != neighbours.second; ++neighbours.first, ++j)
    {
      SP::DynamicalSystem ds2 = (*neighbours.first)->body;
      int ids1 = ds1->number();
      int ids2 = ds2->number();
      int imax = (std::max)(ids1, ids2);
      int imin = (std::min)(ids1, ids2);
      if (ids1 != ids2)
      {
        // is interaction already treated ?
        interPair interpair;
        interpair = std::pair<int, int>(imin, imax);

        if (declaredInteractions.find(interpair)
            == declaredInteractions.end())
        {
          // no, check proximity
          declaredInteractions.insert(interpair);
          ds2->acceptSP(circularFilter);
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

  void visit(SP::SphereLDS ds1)
  {
    // interactions with plans
    for (unsigned int i = 0; i < parent->_plans->size(0); ++i)
    {
      parent->_PlanSphereLDSFilter((*parent->_plans)(i, 0),
                                   (*parent->_plans)(i, 1),
                                   (*parent->_plans)(i, 2),
                                   (*parent->_plans)(i, 3), ds1);
    }

    SP::SiconosVector Q1 = ds1->q();

    double x1 = Q1->getValue(0);
    double y1 = Q1->getValue(1);
    double z1 = Q1->getValue(2);
    SP::Hashed hds1(new Hashed(ds1, (int) floor(x1 / parent->_cellsize),
                               (int) floor(y1 / parent->_cellsize),
                               (int) floor(z1 / parent->_cellsize)));

    // find all other systems that are in the same cells
    std::pair<space_hash::iterator, space_hash::iterator>
    neighbours = parent->_hash_table->equal_range(hds1);

    unsigned int j;
    interPairs declaredInteractions;

    std11::shared_ptr<_SphereLDSFilter> sphereFilter(new _SphereLDSFilter(parent, ds1));

    for (j = 0; neighbours.first != neighbours.second; ++neighbours.first, ++j)
    {
      SP::DynamicalSystem ds2 = (*neighbours.first)->body;
      int ids1 = ds1->number();
      int ids2 = ds2->number();
      int imax = (std::max)(ids1, ids2);
      int imin = (std::min)(ids1, ids2);
      if (ids1 != ids2)
      {
        // is interaction already treated ?
        interPair interpair;
        interpair = std::pair<int, int>(imin, imax);

        if (declaredInteractions.find(interpair)
            == declaredInteractions.end())
        {
          // no, check proximity
          declaredInteractions.insert(interpair);
          ds2->acceptSP(sphereFilter);
        }

      }
    }
  }


  void visit(SP::SphereNEDS ds1)
  {
    // interactions with plans
    for (unsigned int i = 0; i < parent->_plans->size(0); ++i)
    {
      parent->_PlanSphereNEDSFilter((*parent->_plans)(i, 0),
                                    (*parent->_plans)(i, 1),
                                    (*parent->_plans)(i, 2),
                                    (*parent->_plans)(i, 3), ds1);
    }

    SP::SiconosVector Q1 = ds1->q();

    double x1 = Q1->getValue(0);
    double y1 = Q1->getValue(1);
    double z1 = Q1->getValue(2);
    SP::Hashed hds1(new Hashed(ds1, (int) floor(x1 / parent->_cellsize),
                               (int) floor(y1 / parent->_cellsize),
                               (int) floor(z1 / parent->_cellsize)));

    // find all other systems that are in the same cells
    std::pair<space_hash::iterator, space_hash::iterator>
    neighbours = parent->_hash_table->equal_range(hds1);

    unsigned int j;
    interPairs declaredInteractions;

    std11::shared_ptr<_SphereNEDSFilter> sphereFilter(new _SphereNEDSFilter(parent, ds1));

    for (j = 0; neighbours.first != neighbours.second; ++neighbours.first, ++j)
    {
      SP::DynamicalSystem ds2 = (*neighbours.first)->body;
      int ids1 = ds1->number();
      int ids2 = ds2->number();
      int imax = (std::max)(ids1, ids2);
      int imin = (std::min)(ids1, ids2);
      if (ids1 != ids2)
      {
        // is interaction already treated ?
        interPair interpair;
        interpair = std::pair<int, int>(imin, imax);

        if (declaredInteractions.find(interpair)
            == declaredInteractions.end())
        {
          // no, check proximity
          declaredInteractions.insert(interpair);
          ds2->acceptSP(sphereFilter);
        }

      }
    }
  }

  void visit(SP::ExternalBody d)
  {
    d->selfFindInteractions(parent);
  }


};


/* general proximity detection */
void SpaceFilter::buildInteractions(double time)
{
  SP::DynamicalSystemsGraph
  DSG0 = model()->nonSmoothDynamicalSystem()->topology()->dSG(0);

  std11::shared_ptr<_BodyHash>
  hasher(new _BodyHash(*this));

  std11::shared_ptr<_FindInteractions>
  findInteractions(new _FindInteractions(shared_from_this(), time));

  _hash_table->clear();

  // 1: rehash DS
  DynamicalSystemsGraph::VIterator vi, viend;
  for (std11::tie(vi, viend) = DSG0->vertices();
       vi != viend; ++vi)
  {
    // to avoid cast see dual dispatch, visitor pattern
    DSG0->bundle(*vi)->acceptSP(hasher);
  }

  // 2: prox detection
  for (std11::tie(vi, viend) = DSG0->vertices();
       vi != viend; ++vi)
  {
    DSG0->bundle(*vi)->acceptSP(findInteractions);
  }
  model()->simulation()->initOSNS();
}

//std::pair<space_hash::iterator, space_hash::iterator> SpaceFilter::neighbours(SP::Hashed h)
//{
//  return _hash_table.equal_range(h);
//}

bool SpaceFilter::haveNeighbours(SP::Hashed h)
{
  std::pair<space_hash::iterator, space_hash::iterator> neighbours
    = _hash_table->equal_range(h);
  return (neighbours.first != neighbours.second);
}


struct SpaceFilter::_DiskDistance : public SiconosVisitor
{

  using SiconosVisitor::visit;

  double x;
  double y;
  double r;
  double result;

  _DiskDistance(double x, double y, double r)
    : x(x), y(y), r(r)
  {};

  void visit(SP::Disk d)
  {
    double xd = d->q()->getValue(0);
    double yd = d->q()->getValue(1);

    result = (hypot(x - xd, y - yd) - (r + d->getRadius()));
  }
};




/* only for disks at the moment */
double SpaceFilter::minDistance(SP::Hashed h)
{
  std::pair<space_hash::iterator, space_hash::iterator> neighbours
    = _hash_table->equal_range(h);

  SP::SiconosVector q = std11::static_pointer_cast<LagrangianDS>(h->body)->q();

  double dmin = INFINITY;

  {
    SP::Disk disk = std11::static_pointer_cast<Disk>(h->body);

    std11::shared_ptr<_DiskDistance> distance(new _DiskDistance((*q)(0), (*q)(1), disk->getRadius()));

    for (; neighbours.first != neighbours.second; ++neighbours.first)
    {
      (*neighbours.first)->body->acceptSP(distance);

      dmin = (std::min)(dmin, distance->result);
    }
  }

  return dmin;

}

void SpaceFilter::insert(SP::NonSmoothLaw nslaw,
                         long unsigned int id1,
                         long unsigned int id2)
{
  NSLawMatrix& nslaws = *_nslaws;

  // ublas::matrix size type is not the same on 32 bits and 64 bits
  unsigned int id = std::max((unsigned int) id1, (unsigned int) id2);
  nslaws.resize (std::max((unsigned int) nslaws.size1(), (unsigned int) id+1));

  nslaws(id1, id2) = nslaw;
}

SP::NonSmoothLaw SpaceFilter::nslaw(long unsigned int class1, long unsigned class2)
{
  return (*_nslaws)(class1, class2);
}
