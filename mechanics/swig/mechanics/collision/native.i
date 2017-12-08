// -*- c++ -*-
%module(package="siconos.mechanics.collision", directors="1", allprotected="1") native

%include MechanicsBase.i

%fragment("NumPy_Fragments");

// suppress warning
%ignore  STD11::enable_shared_from_this< Hashed >;
%template (sharedHashed) STD11::enable_shared_from_this< Hashed >;

PY_FULL_REGISTER(SpaceFilter);
PY_FULL_REGISTER(SiconosBodies);

// ExternalBody is an astract class and serializers are not generated
// by builder.py
// PY_FULL_REGISTER(ExternalBody);
