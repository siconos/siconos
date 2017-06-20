// -*- c++ -*-
// SWIG interface for Siconos Mechanics/collision
%module(package="collision", directors="1", allprotected="1") base

%include MechanicsBase.i

// Teach SWIG about the SiconosContactorSet base class (std::vector<SiconosContactor>)
class SiconosContactor;
%shared_ptr(std::vector< std11::shared_ptr< SiconosContactor > >);
%template(VectorOfSPSiconosContactor) std::vector< std11::shared_ptr< SiconosContactor > >;

// Ignore some shadowed (redundant for Python) functions
%ignore SiconosShape::setDimensions(SP::SiconosVector dim);

PY_REGISTER_WITHOUT_HEADER(SiconosSphere);
PY_REGISTER_WITHOUT_HEADER(SiconosPlane);
PY_REGISTER_WITHOUT_HEADER(SiconosBox);
PY_REGISTER_WITHOUT_HEADER(SiconosCylinder);
PY_REGISTER_WITHOUT_HEADER(SiconosConvexHull);
PY_REGISTER_WITHOUT_HEADER(SiconosMesh);
PY_REGISTER_WITHOUT_HEADER(SiconosHeightMap);
PY_FULL_REGISTER(SiconosShape);
PY_REGISTER_WITHOUT_HEADER(SiconosContactorSet);
PY_FULL_REGISTER(SiconosContactor);
PY_FULL_REGISTER(BodyDS);
PY_FULL_REGISTER(SiconosCollisionManager);
