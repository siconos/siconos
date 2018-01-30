// -*- c++ -*-
// SWIG interface for Siconos Mechanics/collision
%module(package="siconos.mechanics.collision", directors="1", allprotected="1") base

%include MechanicsBase.i

// Teach SWIG about the SiconosContactorSet base class (std::vector<SiconosContactor>)
class SiconosContactor;
%shared_ptr(std::vector< std11::shared_ptr< SiconosContactor > >);
%template(VectorOfSPSiconosContactor) std::vector< std11::shared_ptr< SiconosContactor > >;

// Result of queries against the collision world is a vector
class SiconosCollisionQueryResult;
%shared_ptr(std::vector< std11::shared_ptr< SiconosCollisionQueryResult > >);
%template(VectorOfSPSiconosCollisionQueryResult) std::vector< std11::shared_ptr<SiconosCollisionQueryResult> >;

// Ignore some shadowed (redundant for Python) functions
%ignore SiconosShape::setDimensions(SP::SiconosVector dim);

PY_REGISTER_WITHOUT_HEADER(SiconosSphere, Mechanics);
PY_REGISTER_WITHOUT_HEADER(SiconosPlane, Mechanics);
PY_REGISTER_WITHOUT_HEADER(SiconosBox, Mechanics);
PY_REGISTER_WITHOUT_HEADER(SiconosCylinder, Mechanics);
PY_REGISTER_WITHOUT_HEADER(SiconosConvexHull, Mechanics);
PY_REGISTER_WITHOUT_HEADER(SiconosMesh, Mechanics);
PY_REGISTER_WITHOUT_HEADER(SiconosHeightMap, Mechanics);
PY_FULL_REGISTER(SiconosShape, Mechanics);
PY_REGISTER_WITHOUT_HEADER(SiconosContactorSet, Mechanics);
PY_FULL_REGISTER(SiconosContactor, Mechanics);
PY_FULL_REGISTER(ContactR, Mechanics);
PY_FULL_REGISTER(BodyDS, Mechanics);
PY_FULL_REGISTER(SiconosCollisionQueryResult, Mechanics);
PY_FULL_REGISTER(SiconosCollisionManager, Mechanics);

%inline
{
  SP::ContactR cast_ContactR(SP::Relation rel)
  {
    return std11::dynamic_pointer_cast<ContactR>(rel);
  };

}
