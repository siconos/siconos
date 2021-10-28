// List of all functions that must be ignored by swig.
// Useful to silent warnings like "Warning 512: Overloaded method ... ignored
//
// See for example http://swig.10945.n7.nabble.com/Python-and-C-operators-just-to-verify-td11222.html
%ignore Siconos::Properties::operator[];
%ignore Siconos::SubProperties::operator[];
%ignore boost::enable_shared_from_this::operator=;
%ignore SiconosGraph::bundle const;
%ignore SiconosGraph::default_color_type const;
%ignore SiconosGraph::color const;
%ignore SiconosGraph::index const;
%ignore SiconosMatrix::operator() const;
%ignore SiconosMatrix::block() const;
%ignore SiconosMatrix::block(unsigned int) const;
%ignore SiconosMatrix::block(unsigned int,unsigned int) const;
%ignore SimpleMatrix::SimpleMatrix(SiconosMatrix const &);
%ignore SimpleMatrix::operator ()(unsigned int,unsigned int) const;
%ignore SimpleMatrix::operator =(SiconosMatrix const &);
%ignore SiconosVector::operator ()(unsigned int) const;
%ignore SiconosVector::operator [](unsigned int) const;
%ignore operator ==(SiconosVector const &,SiconosVector const &);
%ignore BlockVector::begin() const;
%ignore BlockVector::end() const;
%ignore BlockVector::operator ()(unsigned int) const;
%ignore BlockVector::vector(unsigned int) const;
//%ignore BlockVector::operator [](unsigned int) const;
%ignore OneStepIntegrator::dynamicalSystemsBegin() const;
%ignore OneStepIntegrator::dynamicalSystemsEnd() const;
%ignore *::ACCEPT_NONVIRTUAL_VISITORS;

%ignore operator *(double,SiconosVector const &);
%ignore operator *(double,SiconosMatrix const &);
%ignore operator *(SiconosVector const &,double);
%ignore operator *(SiconosMatrix const &,double);
%ignore operator /(SiconosVector const &,double);
%ignore operator /(SiconosMatrix const &,double);
%ignore operator +(SP::SimpleMatrix const,SP::SimpleMatrix const);
%ignore operator +(SiconosMatrix const &,SiconosMatrix const &);
%ignore operator +(SiconosVector const &,SiconosVector const &);
%ignore operator -(SiconosVector const &,SiconosVector const &);
%ignore operator -(SiconosMatrix const &,SiconosMatrix const &);

%ignore xMemory() const;
%ignore normalizeq(SiconosVector &q);

  
