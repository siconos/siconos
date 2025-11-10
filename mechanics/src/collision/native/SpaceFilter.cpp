/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "SpaceFilter.hpp"

#include <boost/container_hash/hash.hpp>  // For hash_combine

#include "Circle.hpp"
#include "CircleCircleR.hpp"
#include "Disk.hpp"
#include "DiskDiskR.hpp"
#include "DiskMovingPlanR.hpp"
#include "DiskPlanR.hpp"
#include "DynamicalSystem.hpp"
#include "DynamicalSystemVisitor.hpp"
#include "ExternalBody.hpp"
#include "Interaction.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "SimulationGraphs.hpp"
#include "SpaceFilter_impl.hpp"
#include "SphereLDS.hpp"
#include "SphereLDSPlanR.hpp"
#include "SphereLDSSphereLDSR.hpp"
#include "SphereNEDS.hpp"
#include "SphereNEDSPlanR.hpp"
#include "SphereNEDSSphereNEDSR.hpp"
#include "Topology.hpp"
// #include<NonSmoothDynamicalSystem.hpp>
// #include <NonSmoothLaw.hpp>

// #include <cmath>
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

// ---- Hashed methods or related functions (declared in SpaceFilter_impl) ----
/* hash is done with encapsulation */
/* needed by boost hash */
bool siconos::collision::native::internal::operator==(std::shared_ptr<Hashed> const& a,
                                                      std::shared_ptr<Hashed> const& b) {
  return (a->i == b->i && a->j == b->j && a->k == b->k);
}

std::size_t siconos::collision::native::internal::hash_value(
    std::shared_ptr<Hashed> const& h) {
  std::size_t seed = 0;
  boost::hash_combine(seed, h->i);
  boost::hash_combine(seed, h->j);
  boost::hash_combine(seed, h->k);
  return seed;
}

siconos::collision::native::SpaceFilter::SpaceFilter()
    : _hash_table(std::make_shared<SpaceHash>()),
      diskdisk_relations(std::make_shared<DiskDiskRDeclaredPool>()),
      diskplan_relations(std::make_shared<DiskPlanRDeclaredPool>()),
      circlecircle_relations(std::make_shared<CircleCircleRDeclaredPool>()) {};

siconos::collision::native::SpaceFilter::SpaceFilter(
    unsigned int bboxfactor, unsigned int cellsize,
    std::shared_ptr<siconos::algebra::SiconosMatrix> plans,
    std::shared_ptr<siconos::collision::native::FMatrix> moving_plans)
    : siconos::simulation::InteractionManager(),
      _bboxfactor(bboxfactor),
      _cellsize(cellsize),
      _plans(plans),
      _moving_plans(moving_plans),
      _hash_table(std::make_shared<SpaceHash>()),
      diskdisk_relations(std::make_shared<DiskDiskRDeclaredPool>()),
      diskplan_relations(std::make_shared<DiskPlanRDeclaredPool>()),
      circlecircle_relations(std::make_shared<CircleCircleRDeclaredPool>()) {};

siconos::collision::native::SpaceFilter::SpaceFilter(
    unsigned int bboxfactor, unsigned int cellsize,
    std::shared_ptr<siconos::algebra::SiconosMatrix> plans)
    : SpaceFilter(bboxfactor, cellsize, plans, nullptr) {}

/* the hashing is done with a visitor */
struct siconos::collision::native::SpaceFilter::BodyHashVisitor_
    : public siconos::modeling::dynamical_systems::Visitor {
 public:
  siconos::simulation::Simulation& sim;
  SpaceFilter& parent;
  BodyHashVisitor_(siconos::simulation::Simulation& s, SpaceFilter& p) : sim(s), parent(p) {};

  using siconos::modeling::dynamical_systems::Visitor::visit;

  void visit(std::shared_ptr<siconos::collision::native::bodies::Disk> pds) override {
    int i, j, imin, imax, jmin, jmax;
    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int)floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int)floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int)floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int)floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i) {
      for (j = jmin; j <= jmax; ++j) {
        parent.insert(pds, i, j, 0);
      }
    }
  };

  void visit(std::shared_ptr<siconos::collision::native::bodies::Circle> pds) override {
    int i, j, imin, imax, jmin, jmax;

    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int)floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int)floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int)floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int)floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i) {
      for (j = jmin; j <= jmax; ++j) {
        parent.insert(pds, i, j, 0);
      }
    }
  }

  void visit(std::shared_ptr<siconos::collision::native::bodies::SphereLDS> pds) override {
    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;

    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int)floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int)floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int)floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int)floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    kmin = (int)floor((pds->getQ(2) - _bboxfactor * pds->getRadius()) / _cellsize);
    kmax = (int)floor((pds->getQ(2) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i) {
      for (j = jmin; j <= jmax; ++j) {
        for (k = kmin; k <= kmax; ++k) {
          parent.insert(pds, i, j, k);
        }
      }
    }
  }

  void visit(std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> pds) override {
    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;

    unsigned int _bboxfactor = parent.bboxfactor();
    unsigned int _cellsize = parent.cellsize();

    imin = (int)floor((pds->getQ(0) - _bboxfactor * pds->getRadius()) / _cellsize);
    imax = (int)floor((pds->getQ(0) + _bboxfactor * pds->getRadius()) / _cellsize);

    jmin = (int)floor((pds->getQ(1) - _bboxfactor * pds->getRadius()) / _cellsize);
    jmax = (int)floor((pds->getQ(1) + _bboxfactor * pds->getRadius()) / _cellsize);

    kmin = (int)floor((pds->getQ(2) - _bboxfactor * pds->getRadius()) / _cellsize);
    kmax = (int)floor((pds->getQ(2) + _bboxfactor * pds->getRadius()) / _cellsize);

    for (i = imin; i <= imax; ++i) {
      for (j = jmin; j <= jmax; ++j) {
        for (k = kmin; k <= kmax; ++k) {
          parent.insert(pds, i, j, k);
        }
      }
    }
  }

  void visit(std::shared_ptr<siconos::collision::native::bodies::ExternalBody> d) override {
    d->selfHash(parent);
  }

  /* ... visit other objects */
};

// --- CircularDSFilterVisitor ---
// A visitor to handle detection for CircularDS
struct siconos::collision::native::SpaceFilter::CircularFilterVisitor_
    : public siconos::modeling::dynamical_systems::Visitor {
  using siconos::modeling::dynamical_systems::Visitor::visit;
  std::shared_ptr<siconos::simulation::Simulation> sim;
  std::shared_ptr<siconos::collision::native::SpaceFilter> parent;
  std::shared_ptr<siconos::collision::native::bodies::CircularDS> ds1;

  CircularFilterVisitor_(std::shared_ptr<siconos::simulation::Simulation> s,
                         std::shared_ptr<siconos::collision::native::SpaceFilter> p,
                         std::shared_ptr<siconos::collision::native::bodies::CircularDS> d1)
      : sim{s}, parent{p}, ds1{d1} {};

  void visit_circular(
      std::shared_ptr<siconos::collision::native::bodies::CircularDS> ds2) const {
    assert(ds1 != ds2);
    auto DSG0 = sim->nonSmoothDynamicalSystem()->topology()->dSG(0);
    assert(DSG0->bundle(DSG0->descriptor(ds1)) == ds1);
    assert(DSG0->bundle(DSG0->descriptor(ds2)) == ds2);
    double r1 = ds1->getRadius();
    double r2 = ds2->getRadius();
    double tol = r1 + r2;
    double x1 = ds1->getQ(0);
    double y1 = ds1->getQ(1);
    double x2 = ds2->getQ(0);
    double y2 = ds2->getQ(1);
    double rmax = std::fmax(r1, r2);
    double rmin = std::fmin(r1, r2);
    double d = std::hypot(x1 - x2, y1 - y2);
    std::shared_ptr<siconos::collision::native::bodies::CircularR> rel = nullptr;
    if (d < rmax) {
      // one inside the other : CircleCircle relation
      if (rmax - (d + rmin) < tol) {
        auto rcandid = parent->circlecircle_relations->find(CircleCircleRDeclared(r1, r2));
        if (rcandid == parent->circlecircle_relations->end()) {
          // a new relation
          rel = std::make_shared<siconos::collision::native::bodies::CircleCircleR>(r1, r2);
          // FIX : this could work with stateless relations.
          // This is not the case: cf LagrangianR.
          //(*(parent->circlecircle_relations))[CircleCircleRDeclared(r1, r2)] = rel;
        } else {
          // get relation from pool
          rel = rcandid->second;
        }
      }
    } else {
      // a DiskDisk relation
      if (d - (r1 + r2) < tol) {
        auto rcandid = parent->diskdisk_relations->find(DiskDiskRDeclared(r1, r2));
        if (rcandid == parent->diskdisk_relations->end()) {
          // a new relation
          rel = std::make_shared<siconos::collision::native::bodies::DiskDiskR>(r1, r2);

          // FIX : this could work with stateless relations.
          // This is not the case cf LagrangianR:
          // parent->diskdisk_relations[DiskDiskRDeclared(r1,r2)] = rel;
        } else {
          // get relation from pool
          rel = rcandid->second;
        }
      }
    }
    if (rel) {
      siconos::simulation::interactions_manager::build_and_link_interaction(ds1, ds2, DSG0,
                                                                            rel, sim, parent);
    } else {
      siconos::simulation::interactions_manager::remove_interaction_if_exists(ds1, ds2, DSG0,
                                                                              sim);
    }
  }

  void visit(std::shared_ptr<siconos::collision::native::bodies::Disk> disk) override {
    visit_circular(disk);
  }

  void visit(std::shared_ptr<siconos::collision::native::bodies::Circle> circle) override {
    visit_circular(circle);
  }

  // do nothing (everything must be done in ExternalBody.findInteractions)
  void visit(std::shared_ptr<siconos::collision::native::bodies::ExternalBody>) override {}
};

/* proximity detection for sphere objects */
struct siconos::collision::native::SpaceFilter::SphereLDSFilterVisitor_
    : public siconos::modeling::dynamical_systems::Visitor {
  using siconos::modeling::dynamical_systems::Visitor::visit;

  std::shared_ptr<siconos::simulation::Simulation> sim;
  std::shared_ptr<siconos::collision::native::SpaceFilter> parent;
  std::shared_ptr<siconos::collision::native::bodies::SphereLDS> ds1;

  SphereLDSFilterVisitor_(std::shared_ptr<siconos::simulation::Simulation> s,
                          std::shared_ptr<siconos::collision::native::SpaceFilter> p,
                          std::shared_ptr<siconos::collision::native::bodies::SphereLDS> ds)
      : sim(s), parent(p), ds1(ds) {};

  void visit(std::shared_ptr<siconos::collision::native::bodies::SphereLDS> ds2) override {
    auto DSG0(sim->nonSmoothDynamicalSystem()->topology()->dSG(0));

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

    if (d < 2 * tol) {
      auto rel =
          std::make_shared<siconos::collision::native::bodies::SphereLDSSphereLDSR>(r1, r2);
      siconos::simulation::interactions_manager::build_and_link_interaction(ds1, ds2, DSG0,
                                                                            rel, sim, parent);
    } else {
      siconos::simulation::interactions_manager::remove_interaction_if_exists(ds1, ds2, DSG0,
                                                                              sim);
    }
  }
};

struct siconos::collision::native::SpaceFilter::SphereNEDSFilterVisitor_
    : public siconos::modeling::dynamical_systems::Visitor {
  using siconos::modeling::dynamical_systems::Visitor::visit;

  std::shared_ptr<siconos::simulation::Simulation> sim;
  std::shared_ptr<siconos::collision::native::SpaceFilter> parent;
  std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> ds1;

  SphereNEDSFilterVisitor_(
      std::shared_ptr<siconos::simulation::Simulation> s,
      std::shared_ptr<siconos::collision::native::SpaceFilter> p,
      std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> neds)
      : sim(s), parent(p), ds1(neds) {};

  void visit(std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> ds2) override {
    auto DSG0 = sim->nonSmoothDynamicalSystem()->topology()->dSG(0);

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

    if (d < 2 * tol) {
      auto rel =
          std::make_shared<siconos::collision::native::bodies::SphereNEDSSphereNEDSR>(r1, r2);
      siconos::simulation::interactions_manager::build_and_link_interaction(ds1, ds2, DSG0,
                                                                            rel, sim, parent);
    } else {
      // is interaction in graph ?
      bool found = false;
      siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
      for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds1)); oei != oeiend;
           ++oei) {
        if (DSG0->bundle(DSG0->target(*oei)) == ds2) {
          found = true;
          break;
        }
      }

      if (found) {
        sim->unlink(DSG0->bundle(*oei));
      }
    }
  }
};

/* disk plan relation comparison */
struct siconos::collision::native::SpaceFilter::IsSameDiskPlanRVisitor_
    : public siconos::modeling::relations::Visitor {
  using siconos::modeling::relations::Visitor::visit;

  std::shared_ptr<siconos::simulation::Simulation> sim;
  std::shared_ptr<siconos::collision::native::SpaceFilter> parent;
  double A, B, C, r, xCenter, yCenter, width;
  bool flag;
  IsSameDiskPlanRVisitor_(std::shared_ptr<siconos::simulation::Simulation> s,
                          std::shared_ptr<siconos::collision::native::SpaceFilter> p, double A,
                          double B, double C, double r, double xCenter, double yCenter,
                          double width)
      : sim(s),
        parent(p),
        A(A),
        B(B),
        C(C),
        r(r),
        xCenter(xCenter),
        yCenter(yCenter),
        width(width),
        flag(false) {};

  void visit(const siconos::collision::native::bodies::DiskDiskR&) override { flag = false; };

  void visit(const siconos::collision::native::bodies::CircleCircleR&) override {
    flag = false;
  };

  void visit(const siconos::collision::native::bodies::DiskMovingPlanR&) override {
    flag = false;
  };

  void visit(const siconos::modeling::LagrangianScleronomousR&) override { flag = false; };

  void visit(const siconos::collision::native::bodies::DiskPlanR& rel) override {
    flag = rel.equal(A, B, C, r, xCenter, yCenter, width);
  };
};

struct siconos::collision::native::SpaceFilter::IsSameDiskMovingPlanRVisitor_
    : public siconos::modeling::relations::Visitor {
  using siconos::modeling::relations::Visitor::visit;

  std::shared_ptr<siconos::collision::native::SpaceFilter> parent;
  siconos::modeling::func_prototypes::FunctionS_S AF, BF, CF;
  double r;
  bool flag;
  IsSameDiskMovingPlanRVisitor_(std::shared_ptr<siconos::collision::native::SpaceFilter> p,
                                const siconos::modeling::func_prototypes::FunctionS_S& AF,
                                const siconos::modeling::func_prototypes::FunctionS_S& BF,
                                const siconos::modeling::func_prototypes::FunctionS_S& CF,
                                double r)
      : parent(p), AF(AF), BF(BF), CF(CF), r(r), flag(false) {};

  void visit(const siconos::collision::native::bodies::DiskDiskR&) override { flag = false; };

  void visit(const siconos::collision::native::bodies::CircleCircleR&) override {
    flag = false;
  };

  void visit(const siconos::collision::native::bodies::DiskPlanR&) override { flag = false; };

  void visit(const siconos::modeling::LagrangianScleronomousR&) override { flag = false; }

  void visit(const siconos::collision::native::bodies::DiskMovingPlanR& rel) override {
    flag = rel.equal(AF, BF, CF, r);
  };
};

/* sphere plan relation comparison */
struct siconos::collision::native::SpaceFilter::IsSameSpherePlanRVisitor_
    : public siconos::modeling::relations::Visitor {
  using siconos::modeling::relations::Visitor::visit;

  std::shared_ptr<siconos::collision::native::SpaceFilter> parent;
  double A, B, C, D, r;
  bool flag;
  IsSameSpherePlanRVisitor_(std::shared_ptr<siconos::collision::native::SpaceFilter> p,
                            double A, double B, double C, double D, double r)
      : parent(p), A(A), B(B), C(C), D(D), r(r), flag(false) {};

  void visit(const siconos::collision::native::bodies::SphereLDSSphereLDSR&) override {
    flag = false;
  };

  void visit(const siconos::collision::native::bodies::SphereNEDSSphereNEDSR&) override {
    flag = false;
  };

  void visit(const siconos::collision::native::bodies::SphereLDSPlanR& rel) override {
    flag = rel.equal(A, B, C, D, r);
  };

  void visit(const siconos::collision::native::bodies::SphereNEDSPlanR& rel) override {
    flag = rel.equal(A, B, C, D, r);
  };
};

/* proximity detection between circular object and plans */
void siconos::collision::native::SpaceFilter::_PlanCircularFilter(
    std::shared_ptr<siconos::simulation::Simulation> sim, double A, double B, double C,
    double xCenter, double yCenter, double width,
    std::shared_ptr<siconos::collision::native::bodies::CircularDS> ds) {
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  auto DSG0 = sim->nonSmoothDynamicalSystem()->topology()->dSG(0);

  IsSameDiskPlanRVisitor_ isSameDiskPlanR{
      sim,     std::dynamic_pointer_cast<SpaceFilter>(shared_from_this()),
      A,       B,
      C,       r,
      xCenter, yCenter,
      width};

  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  auto relp = std::make_shared<siconos::collision::native::bodies::DiskPlanR>(
      r, A, B, C, xCenter, yCenter, width);

  if (relp->distance(ds->getQ(0), ds->getQ(1), ds->getRadius()) < tol) {
    // is interaction in graph ?
    bool found = false;
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameDiskPlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameDiskPlanR.flag) {
        found = true;
        break;
      }
    }

    if (!found)
    // no
    {
      auto nslaw = nonSmoothLaw(DSG0->groupId[DSG0->descriptor(ds)],
                                DSG0->groupId[DSG0->descriptor(ds)]);

      auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relp);
      DEBUG_PRINTF("insert interaction : %d\n", inter->number());
      sim->link(inter, ds);
    }
  } else {
    // is interaction in graph ?
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameDiskPlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameDiskPlanR.flag) {
        DEBUG_PRINTF("remove interaction : %d\n", DSG0->bundle(*oei)->number());

        sim->unlink(DSG0->bundle(*oei));
        break;
      }
    }
  }
}

/* proximity detection between circular object and plans */
void siconos::collision::native::SpaceFilter::_MovingPlanCircularFilter(
    std::shared_ptr<siconos::simulation::Simulation> sim, unsigned int i,
    std::shared_ptr<siconos::collision::native::bodies::CircularDS> ds, double time) {
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  auto DSG0 = sim->nonSmoothDynamicalSystem()->topology()->dSG(0);

  IsSameDiskMovingPlanRVisitor_ isSameDiskMovingPlanR{
      std::dynamic_pointer_cast<SpaceFilter>(shared_from_this()), (*_moving_plans)(i, 0),
      (*_moving_plans)(i, 1), (*_moving_plans)(i, 2), r};

  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  auto relp = std::make_shared<siconos::collision::native::bodies::DiskMovingPlanR>(r);

  relp->setComputeAFunction((*_moving_plans)(i, 0));
  relp->setComputeBFunction((*_moving_plans)(i, 1));
  relp->setComputeCFunction((*_moving_plans)(i, 2));
  relp->setComputeAdotFunction((*_moving_plans)(i, 3));
  relp->setComputeBdotFunction((*_moving_plans)(i, 4));
  relp->setComputeCdotFunction((*_moving_plans)(i, 5));

  relp->init(time);

  if (relp->distance(ds->getQ(0), ds->getQ(1), ds->getRadius()) < tol)

  {
    // is interaction in graph ?
    bool found = false;
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameDiskMovingPlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameDiskMovingPlanR.flag) {
        found = true;
        break;
      }
    }
    if (!found)
    // no
    {
      auto nslaw = nonSmoothLaw(DSG0->groupId[DSG0->descriptor(ds)],
                                DSG0->groupId[DSG0->descriptor(ds)]);

      auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relp);
      sim->link(inter, ds);
    }
  } else {
    // is interaction in graph ?
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameDiskMovingPlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameDiskMovingPlanR.flag) {
        sim->unlink(DSG0->bundle(*oei));
        break;
      }
    }
  }
}

/* proximity detection between circular object and plans */
void siconos::collision::native::SpaceFilter::_PlanSphereLDSFilter(
    std::shared_ptr<siconos::simulation::Simulation> sim, double A, double B, double C,
    double D, std::shared_ptr<siconos::collision::native::bodies::SphereLDS> ds) {
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  auto DSG0 = sim->nonSmoothDynamicalSystem()->topology()->dSG(0);

  IsSameSpherePlanRVisitor_ isSameSpherePlanR{
      std::dynamic_pointer_cast<SpaceFilter>(shared_from_this()), A, B, C, D, r};

  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  auto relp =
      std::make_shared<siconos::collision::native::bodies::SphereLDSPlanR>(r, A, B, C, D);
  if (relp->distance(ds->getQ(0), ds->getQ(1), ds->getQ(2), ds->getRadius()) < tol)

  {
    // is interaction in graph ?
    bool found = false;
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameSpherePlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameSpherePlanR.flag) {
        found = true;
        break;
      }
    }
    if (!found)
    // no
    {
      auto nslaw = nonSmoothLaw(DSG0->groupId[DSG0->descriptor(ds)],
                                DSG0->groupId[DSG0->descriptor(ds)]);

      auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relp);
      sim->link(inter, ds);
    }
  } else {
    // is interaction in graph ?
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameSpherePlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameSpherePlanR.flag) {
        sim->unlink(DSG0->bundle(*oei));
        break;
      }
    }
  }
}

// note : all PlanObject should be merged
void siconos::collision::native::SpaceFilter::_PlanSphereNEDSFilter(
    std::shared_ptr<siconos::simulation::Simulation> sim, double A, double B, double C,
    double D, std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> ds) {
  double r = ds->getRadius();

  /* tolerance */
  double tol = r;

  auto DSG0 = sim->nonSmoothDynamicalSystem()->topology()->dSG(0);

  IsSameSpherePlanRVisitor_ isSameSpherePlanR{
      std::dynamic_pointer_cast<SpaceFilter>(shared_from_this()), A, B, C, D, r};

  // all DS must be in DS graph
  assert(DSG0->bundle(DSG0->descriptor(ds)) == ds);
  auto relp =
      std::make_shared<siconos::collision::native::bodies::SphereNEDSPlanR>(r, A, B, C, D);
  if (relp->distance(ds->getQ(0), ds->getQ(1), ds->getQ(2), ds->getRadius()) < tol)

  {
    // is interaction in graph ?
    bool found = false;
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameSpherePlanR);
      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameSpherePlanR.flag) {
        found = true;
        break;
      }
    }
    if (!found)
    // no
    {
      auto nslaw = nonSmoothLaw(DSG0->groupId[DSG0->descriptor(ds)],
                                DSG0->groupId[DSG0->descriptor(ds)]);

      auto inter = std::make_shared<siconos::modeling::Interaction>(nslaw, relp);
      sim->link(inter, ds);
    }
  } else {
    // is interaction in graph ?
    siconos::graphs::DynamicalSystemsGraph::OEIterator oei, oeiend;
    for (std::tie(oei, oeiend) = DSG0->out_edges(DSG0->descriptor(ds)); oei != oeiend; ++oei) {
      DSG0->bundle(*oei)->relation()->accept(isSameSpherePlanR);

      if (DSG0->bundle(DSG0->target(*oei)) == ds && isSameSpherePlanR.flag) {
        sim->unlink(DSG0->bundle(*oei));
        break;
      }
    }
  }
}

/* insertion */
void siconos::collision::native::SpaceFilter::insert(
    std::shared_ptr<siconos::collision::native::bodies::Disk> ds, int i, int j, int k) {
  auto hashed = std::make_shared<siconos::collision::native::internal::Hashed>(
      std::static_pointer_cast<siconos::modeling::DynamicalSystem>(ds), i, j);
  _hash_table->insert(hashed);
}

void siconos::collision::native::SpaceFilter::insert(
    std::shared_ptr<siconos::collision::native::bodies::Circle> ds, int i, int j, int k) {
  auto hashed = std::make_shared<siconos::collision::native::internal::Hashed>(
      std::static_pointer_cast<siconos::modeling::DynamicalSystem>(ds), i, j);
  _hash_table->insert(hashed);
}

void siconos::collision::native::SpaceFilter::insert(
    std::shared_ptr<siconos::collision::native::bodies::SphereLDS> ds, int i, int j, int k) {
  auto hashed = std::make_shared<siconos::collision::native::internal::Hashed>(ds, i, j, k);
  _hash_table->insert(hashed);
}

void siconos::collision::native::SpaceFilter::insert(
    std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> ds, int i, int j, int k) {
  auto hashed = std::make_shared<siconos::collision::native::internal::Hashed>(ds, i, j, k);
  _hash_table->insert(hashed);
}

void siconos::collision::native::SpaceFilter::insert(
    std::shared_ptr<siconos::collision::native::internal::Hashed> hashed) {
  _hash_table->insert(hashed);
}

/* insert other objects */

/* dynamical systems proximity detection */

bool siconos::collision::native::operator==(
    siconos::collision::native::DiskPlanRDeclared const& a,
    siconos::collision::native::DiskPlanRDeclared const& b) {
  return ((a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3] && a[4] == b[4] &&
           a[5] == b[5]));
}

struct siconos::collision::native::SpaceFilter::FindInteractionsVisitor_
    : public siconos::modeling::dynamical_systems::Visitor {
  using siconos::modeling::dynamical_systems::Visitor::visit;

  using InterPairs =
      std::unordered_multiset<std::pair<int, int>, boost::hash<std::pair<int, int>>>;

  std::shared_ptr<siconos::simulation::Simulation> sim{nullptr};
  std::shared_ptr<siconos::collision::native::SpaceFilter> parent{nullptr};
  double time{0.};

  FindInteractionsVisitor_(std::shared_ptr<siconos::simulation::Simulation> s,
                           std::shared_ptr<siconos::collision::native::SpaceFilter> p,
                           double time)
      : sim(s), parent(p), time(time) {};

  void visit_circular(std::shared_ptr<siconos::collision::native::bodies::CircularDS> ds1) {
    assert(parent->_plans->rows() > 0);

    // interactions with plans

    if (parent->_plans) {
      siconos::algebra::print(*parent->_plans);
      for (siconos::algebra::Index i = 0; i < parent->_plans->rows(); ++i) {
        parent->_PlanCircularFilter(sim, (*parent->_plans)(i, 0), (*parent->_plans)(i, 1),
                                    (*parent->_plans)(i, 2), (*parent->_plans)(i, 3),
                                    (*parent->_plans)(i, 4), (*parent->_plans)(i, 5), ds1);
      }
    }

    if (parent->_moving_plans) {
      for (unsigned int i = 0; i < parent->_moving_plans->size1(); ++i) {
        parent->_MovingPlanCircularFilter(sim, i, ds1, time);
      }
    }

    auto Q1 = ds1->q();

    auto x1 = (*Q1)(0);
    auto y1 = (*Q1)(1);
    auto hds1 = std::make_shared<siconos::collision::native::internal::Hashed>(
        std::static_pointer_cast<siconos::modeling::DynamicalSystem>(ds1),
        (int)floor(x1 / parent->_cellsize), (int)floor(y1 / parent->_cellsize));

    // find all other systems that are in the same cells
    auto neighbours = parent->_hash_table->equal_range(hds1);

    InterPairs declaredInteractions;
    CircularFilterVisitor_ circularFilter(sim, parent, ds1);
    for (unsigned int j = 0; neighbours.first != neighbours.second; ++neighbours.first, ++j) {
      auto ds2 = (*(neighbours.first))->body;
      auto ids1 = ds1->number();
      auto ids2 = ds2->number();
      auto imax = (std::max)(ids1, ids2);
      auto imin = (std::min)(ids1, ids2);
      if (ids1 != ids2) {
        // is interaction already treated ?
        auto interpair = std::pair<int, int>(imin, imax);

        if (declaredInteractions.find(interpair) == declaredInteractions.end()) {
          // no, check proximity
          declaredInteractions.insert(interpair);
          ds2->acceptSP(circularFilter);
        }
      }
    }
  };

  void visit(std::shared_ptr<siconos::collision::native::bodies::Circle> circle) override {
    visit_circular(circle);
  };

  void visit(std::shared_ptr<siconos::collision::native::bodies::Disk> disk) override {
    visit_circular(disk);
  };

  void visit(std::shared_ptr<siconos::collision::native::bodies::SphereLDS> ds1) override {
    // interactions with plans
    for (siconos::algebra::Index i = 0; i < parent->_plans->rows(); ++i) {
      parent->_PlanSphereLDSFilter(sim, (*parent->_plans)(i, 0), (*parent->_plans)(i, 1),
                                   (*parent->_plans)(i, 2), (*parent->_plans)(i, 3), ds1);
    }

    auto Q1 = ds1->q();

    double x1 = (*Q1)(0);
    double y1 = (*Q1)(1);
    double z1 = (*Q1)(2);
    auto hds1 = std::make_shared<siconos::collision::native::internal::Hashed>(
        ds1, (int)floor(x1 / parent->_cellsize), (int)floor(y1 / parent->_cellsize),
        (int)floor(z1 / parent->_cellsize));

    // find all other systems that are in the same cells
    auto neighbours = parent->_hash_table->equal_range(hds1);

    InterPairs declaredInteractions;

    SphereLDSFilterVisitor_ sphereFilter{sim, parent, ds1};
    for (unsigned int j = 0; neighbours.first != neighbours.second; ++neighbours.first, ++j) {
      auto ds2 = (*neighbours.first)->body;
      int ids1 = ds1->number();
      int ids2 = ds2->number();
      int imax = (std::max)(ids1, ids2);
      int imin = (std::min)(ids1, ids2);
      if (ids1 != ids2) {
        // is interaction already treated ?
        auto interpair = std::pair<int, int>(imin, imax);

        if (declaredInteractions.find(interpair) == declaredInteractions.end()) {
          // no, check proximity
          declaredInteractions.insert(interpair);
          ds2->acceptSP(sphereFilter);
        }
      }
    }
  }

  void visit(std::shared_ptr<siconos::collision::native::bodies::SphereNEDS> ds1) override {
    // interactions with plans
    for (siconos::algebra::Index i = 0; i < parent->_plans->rows(); ++i) {
      parent->_PlanSphereNEDSFilter(sim, (*parent->_plans)(i, 0), (*parent->_plans)(i, 1),
                                    (*parent->_plans)(i, 2), (*parent->_plans)(i, 3), ds1);
    }

    auto Q1 = ds1->q();

    double x1 = (*Q1)(0);
    double y1 = (*Q1)(1);
    double z1 = (*Q1)(2);
    auto hds1 = std::make_shared<siconos::collision::native::internal::Hashed>(
        ds1, (int)floor(x1 / parent->_cellsize), (int)floor(y1 / parent->_cellsize),
        (int)floor(z1 / parent->_cellsize));

    // find all other systems that are in the same cells
    auto neighbours = parent->_hash_table->equal_range(hds1);

    InterPairs declaredInteractions;

    SphereNEDSFilterVisitor_ sphereFilter{sim, parent, ds1};
    for (unsigned int j = 0; neighbours.first != neighbours.second; ++neighbours.first, ++j) {
      auto ds2 = (*neighbours.first)->body;
      int ids1 = ds1->number();
      int ids2 = ds2->number();
      int imax = (std::max)(ids1, ids2);
      int imin = (std::min)(ids1, ids2);
      if (ids1 != ids2) {
        // is interaction already treated ?
        auto interpair = std::pair<int, int>(imin, imax);

        if (declaredInteractions.find(interpair) == declaredInteractions.end()) {
          // no, check proximity
          declaredInteractions.insert(interpair);
          ds2->acceptSP(sphereFilter);
        }
      }
    }
  }

  void visit(std::shared_ptr<siconos::collision::native::bodies::ExternalBody> d) override {
    d->selfFindInteractions(parent);
  }
};

/* general proximity detection */
void siconos::collision::native::SpaceFilter::updateInteractions(
    std::shared_ptr<siconos::simulation::Simulation> sim) {
  double time = sim->getTk();
  auto DSG0 = sim->nonSmoothDynamicalSystem()->topology()->dSG(0);

  BodyHashVisitor_ hasher{*sim, *this};
  FindInteractionsVisitor_ findInteractions{
      sim, std::dynamic_pointer_cast<SpaceFilter>(shared_from_this()), time};

  _hash_table->clear();

  // 1: rehash DS
  siconos::graphs::DynamicalSystemsGraph::VIterator vi, viend;
  for (std::tie(vi, viend) = DSG0->vertices(); vi != viend; ++vi) {
    // to avoid cast see dual dispatch, visitor pattern
    DSG0->bundle(*vi)->acceptSP(hasher);

    // --> this will fill the hash_table with hashed dynamical systems through accept/visit
    // calls.
  }

  // 2: prox detection
  for (std::tie(vi, viend) = DSG0->vertices(); vi != viend; ++vi) {
    DSG0->bundle(*vi)->acceptSP(findInteractions);
  }
  // model()->simulation()->initializeOneStepNSProblem();
}

bool siconos::collision::native::SpaceFilter::haveNeighbours(
    std::shared_ptr<siconos::collision::native::internal::Hashed> h) {
  auto neighbours = _hash_table->equal_range(h);
  return (neighbours.first != neighbours.second);
}

struct siconos::collision::native::SpaceFilter::DiskDistanceVisitor_
    : public siconos::modeling::dynamical_systems::Visitor {
  using siconos::modeling::dynamical_systems::Visitor::visit;

  double x;
  double y;
  double r;
  double result;

  DiskDistanceVisitor_(double x, double y, double r) : x(x), y(y), r(r) {};

  void visit(std::shared_ptr<siconos::collision::native::bodies::Disk> d) override {
    double xd = (*d->q())(0);
    double yd = (*d->q())(1);

    result = (hypot(x - xd, y - yd) - (r + d->getRadius()));
  }
};

/* only for disks at the moment */
double siconos::collision::native::SpaceFilter::minDistance(
    std::shared_ptr<siconos::collision::native::internal::Hashed> h) {
  auto neighbours = _hash_table->equal_range(h);

  auto q = std::static_pointer_cast<siconos::modeling::LagrangianDS>(h->body)->q();

  double dmin = std::numeric_limits<double>::infinity();

  {
    auto disk = std::static_pointer_cast<siconos::collision::native::bodies::Disk>(h->body);

    DiskDistanceVisitor_ distance{(*q)(0), (*q)(1), disk->getRadius()};

    for (; neighbours.first != neighbours.second; ++neighbours.first) {
      (*neighbours.first)->body->acceptSP(distance);

      dmin = (std::min)(dmin, distance.result);
    }
  }

  return dmin;
}

void siconos::collision::native::SpaceFilter::insertLine(double a, double b, double c) {
  size_t row;
  if (!_plans) {
    _plans = std::make_shared<siconos::algebra::SiconosMatrix>(1, 6);
    _plans->setZero();
    row = 0;
  } else {
    _plans->conservativeResize(_plans->rows() + 1, 6);
    row = _plans->rows() - 1;
  }
  (*_plans)(row, 0) = a;
  (*_plans)(row, 1) = b;
  (*_plans)(row, 2) = c;
}
