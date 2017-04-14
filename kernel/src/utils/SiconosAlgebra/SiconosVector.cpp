/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

#include "SiconosConfig.h"

#include "SiconosAlgebraTypeDef.hpp"

#include <boost/numeric/ublas/io.hpp>            // for >>
//#include <boost/numeric/ublas/vector_proxy.hpp>  // for project
#include <boost/numeric/ublas/vector_sparse.hpp>


#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace bindings = boost::numeric::bindings::blas;

#include "SiconosVectorStorage.hpp"

#include "SimpleMatrix.hpp"

#include "ioVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosAlgebra.hpp"
#include "BlockVector.hpp"
//#define DEBUG_MESSAGES
#include "debug.h"
#include "TypeName.hpp"

#include <boost/move/unique_ptr.hpp>

#include <boost/tuple/tuple.hpp>

#undef VISITOR_CLASSES
#define VISITOR_CLASSES() \
  REGISTER(DenseVectStorage) \
  REGISTER(SparseVectStorage)

#include <VisitorMaker.hpp>

#include <iostream>
#include <typeinfo>

#include <cmath>        // std::exp(double)
#include <algorithm>    // std::transform


using namespace Experimental;

/* a visitor with some const arguments */
template<typename P1, typename P2=empty, typename P3=empty, typename P4=empty>
struct ParamVisitor : public SiconosVisitor
{
  typedef boost::tuple<const P1&, const P2&, const P3&, const P4&> arguments_type;
  const arguments_type& _param;

  ParamVisitor(const arguments_type& args) : _param(args) {};

  virtual ~ParamVisitor() {};

  const P1& param1() const { return boost::get<0>(_param); };
  const P2& param2() const { return boost::get<1>(_param); };
  const P3& param3() const { return boost::get<2>(_param); };
  const P4& param4() const { return boost::get<3>(_param); };
};


template<>
struct ParamVisitor<empty, empty, empty, empty> : public SiconosVisitor
{
  virtual ~ParamVisitor() {};
};

template<typename P>
struct ParamVisitor<P, empty, empty, empty> : public SiconosVisitor
{
  typedef const P arguments_type;
  const arguments_type& _param;

  ParamVisitor(const arguments_type& a) : _param(a) {};

  virtual ~ParamVisitor() {};

  const P& param() const { return _param; };
  const P& param1() const { return _param; };
};

template<typename P1, typename P2>
struct ParamVisitor<P1, P2, empty, empty> : public SiconosVisitor
{
  typedef boost::tuple<const P1&, const P2&> arguments_type;
  const arguments_type& _param;

  ParamVisitor(const arguments_type& args) : _param(args) {};

  virtual ~ParamVisitor() {};

  const P1& param1() const { return boost::get<0>(_param); };
  const P2& param2() const { return boost::get<1>(_param); };

};

template<typename P1, typename P2, typename P3>
struct ParamVisitor<P1, P2, P3, empty> : public SiconosVisitor
{
  typedef boost::tuple<const P1&, const P2&, const P3&> arguments_type;
  const arguments_type& _param;

  ParamVisitor(const arguments_type& args) : _param(args) {};

  virtual ~ParamVisitor() {};

  const P1& param1() const { return boost::get<0>(_param); };
  const P2& param2() const { return boost::get<1>(_param); };
  const P3& param3() const { return boost::get<2>(_param); };

};

template<typename C, typename P1, typename P2=empty, typename P3=empty>
struct ParamModifier : public SiconosVisitor
{
  const typedef boost::tuple<const P1&, const P2&, const P3&, C&> arguments_type;
  const arguments_type& _param;

  ParamModifier(arguments_type& args) : _param(args) {};

  virtual ~ParamModifier() {};

  C& modifiable() { return boost::get<3>(_param); };
  const P1& param1() const { return boost::get<0>(_param); };
  const P2& param2() const { return boost::get<1>(_param); };
  const P3& param3() const { return boost::get<2>(_param); };
};

template<typename C, typename P1, typename P2>
struct ParamModifier<C, P1, P2, empty> : public SiconosVisitor
{
  typedef boost::tuple<const P1&, const P2&, C&> arguments_type;
  const arguments_type& _param;

  ParamModifier(const arguments_type& args) : _param(args) {};

  virtual ~ParamModifier() {};

  C& modifiable() { return boost::get<2>(_param); };
  const P1& param1() const { return boost::get<0>(_param); };
  const P2& param2() const { return boost::get<1>(_param); };
};

template<typename C, typename P1>
struct ParamModifier<C, P1, empty, empty> : public SiconosVisitor
{
  typedef boost::tuple<const P1&, C&> arguments_type;
  const arguments_type& _param;

  ParamModifier(const arguments_type& args) : _param(args) {};

  virtual ~ParamModifier() {};

  C& modifiable() { return boost::get<1>(_param); };
  const P1& param1() const { return boost::get<0>(_param); };
};



template<typename C>
struct SimpleModifier : public SiconosVisitor
{
  typedef C arguments_type;
  arguments_type& _param;

  SimpleModifier(arguments_type& object) : _param(object) {};
  arguments_type& modifiable() { return _param; };

  virtual ~SimpleModifier() {};
};

#undef REGISTER
#define REGISTER(X) X,
template <typename V, typename P=empty>
struct make_visitor : public Visitor<Classes< VISITOR_CLASSES() empty>, V>::Make
{

  typedef typename Visitor<Classes< VISITOR_CLASSES() empty>, V>::Make base_type;

  make_visitor()  {};

  make_visitor(const typename base_type::arguments_type& args) :
    base_type(args) {};

};

template<typename V, typename T>
void apply_visitor(T& obj)
{
  make_visitor<V> visitor = make_visitor<V>();
  obj.accept(visitor);
};

template<typename V, typename T, typename P>
void apply_visitor(T& obj, const P& param)
{
  const typename make_visitor<V, P>::arguments_type& args = typename make_visitor<V, P>::arguments_type(param);
  make_visitor<V, P> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename T>
void apply_visitor(T& obj, unsigned int param)
{
  make_visitor<V, unsigned int> visitor(param);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1, typename P2>
void apply_visitor(T& obj, const P1& param1, const P2& param2)
{
  typedef typename V::arguments_type P;
  const typename make_visitor<V, P>::arguments_type& args = typename make_visitor<V, P>::arguments_type(param1, param2);
  make_visitor<V, P> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1, typename P2, typename P3>
void apply_visitor(T& obj, const P1& param1, const P2& param2, const P3& param3)
{
  typedef typename V::arguments_type P;
  typename make_visitor<V, P>::arguments_type& args = typename make_visitor<V, P>::arguments_type(param1, param2, param3);
  make_visitor<V, P> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1, typename P2, typename P3, typename P4>
void apply_visitor(T& obj, const P1& param1, const P2& param2, const P3& param3, const P4& param4)
{
  typedef typename V::arguments_type P;
  typename make_visitor<V, P>::arguments_type& args = typename make_visitor<V, P>::arguments_type(param1, param2, param3, param4);
  make_visitor<V, P> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename R, typename T>
R apply_visitor(T& obj)
{
  make_visitor<V> visitor;
  obj.accept(visitor);
  return visitor.answer;
};

template<typename V, typename R, typename P, typename T>
R apply_visitor(T& obj, const P& param)
{
  make_visitor<V, P> visitor(param);
  obj.accept(visitor);
  return visitor.answer;
};

template<typename V, typename R, typename T>
R apply_visitor(T& obj, unsigned int param)
{
  make_visitor<V, unsigned int> visitor(param);
  obj.accept(visitor);
  return visitor.answer;
};


/* modifiers equivalents */
template <typename V>
struct make_modifier : public Modifier<Classes< VISITOR_CLASSES() empty>, V>::Make
{

  typedef typename Modifier<Classes< VISITOR_CLASSES() empty>, V>::Make base_type;
  typedef typename base_type::arguments_type arguments_type;

  make_modifier(arguments_type& p) :
    base_type(p) {};

  make_modifier() {};
};


template<typename V, typename P>
typename Visitor<Classes< VISITOR_CLASSES() empty>, V>::Make xmake_modifier(typename Visitor<Classes< VISITOR_CLASSES() empty>, V>::Make::arguments_type& p) { return typename Visitor<Classes< VISITOR_CLASSES() empty>, V>::Make(p);};


template<typename V, typename T>
void apply_modifier(T& obj)
{
  make_modifier<V> modifier = make_modifier<V>();
  obj.accept_modifier(modifier);
};

template<typename V, typename T, typename P>
void apply_modifier(T& obj, const P& param)
{
  make_modifier<V> visitor = make_modifier<V>(param);
  obj.accept_modifier(visitor);
};

template<typename V, typename T, typename P>
void apply_param_modifier(const T& obj, P& param)
{
  make_modifier<V> modifier = make_modifier<V>(param);
  obj.accept(modifier);
};

template<typename V, typename T, typename P1, typename P2>
void apply_modifier(T& obj, const P1& param1, P2& param2)
{
  typename make_modifier<V>::arguments_type args = typename make_modifier<V>::arguments_type(param1, param2);
  make_modifier<V> modifier(args);
  obj.accept(modifier);
};

template<typename V, typename T, typename P1, typename P2, typename P3>
void apply_modifier(T& obj, const P1& param1, const P2& param2, P3& param3)
{
  typename make_modifier<V>::arguments_type args = typename make_modifier<V>::arguments_type(param1, param2, param3);
  make_modifier<V> modifier(args);
  obj.accept(modifier);
};

template<typename V, typename T, typename P1, typename P2, typename P3, typename P4>
void apply_modifier(const T& obj, const P1& param1, const P2& param2, const P3& param3, P4& param4)
{
  typename make_modifier<V>::arguments_type args = typename make_modifier<V>::arguments_type(param1, param2, param3, param4);
  make_modifier<V> modifier(args);
  obj.accept(modifier);
};

template<typename V, typename R, typename T>
R apply_modifier(T& obj)
{
  make_modifier<V> modifier;
  obj.accept_modifier(modifier);
  return modifier.answer;
};

template<typename V, typename R, typename P, typename T>
R apply_modifier(T& obj, const P& param)
{
  make_modifier<V> modifier(param);
  obj.accept_modifier(modifier);
  return modifier.answer;
};
/* --- */


struct Resize : public ParamVisitor<size_t>
{

  Resize(size_t size) : ParamVisitor(size) {};

  template<typename T>
  void operator() (T& v)
  {
    return v.internal_data.resize(this->param());
  };
};

struct Fill : public ParamVisitor<const double>
{

  Fill(const double& value) : ParamVisitor(value) {};

  template<typename T>
  void operator() (T& v)
  {
    for (unsigned int i = 0; i < v.internal_data.size(); ++i)
    {
      v.internal_data[i] = this->param();
    }
  }
};

template<>
void Fill::operator()<DenseVect> (DenseVect& v)
  {
    bindings::set(this->param(), v);
  }

struct Size : public SiconosVisitor
{
  size_t answer;
  template<typename T>
  void operator() (const T& v)
  {
    answer = v.internal_data.size();
  };
};

struct Zero : public SiconosVisitor
{
  template<typename T>
  void operator () (T& storage)
  {
    storage.internal_data *= 0.;
  }
};

template<>
void Zero::operator()<DenseVectStorage> (DenseVectStorage& dv)
{
  bindings::scal(0.0, dv.internal_data);
};

struct NormInf : public SiconosVisitor
{
  double answer;

  template<typename T>
  void operator () (const T& storage)
  {
    answer = norm_inf(storage.internal_data);
  };
};

struct Norm2 : public SiconosVisitor
{
  double answer;

  template<typename T>
  void operator () (const T& storage)
  {
    answer = ublas::norm_2(storage.internal_data);
  };
};

struct Sum : public SiconosVisitor
{
  double answer;

  template<typename T>
  void operator () (const T& storage)
  {
    answer = sum(storage.internal_data);
  };
};

struct Display : public ParamVisitor<unsigned int>
{

  Display(unsigned int precision) : ParamVisitor(precision) {};

  template<typename T>
  void operator () (const T& storage)
  {
    std::cout.setf(std::ios::scientific);
    std::cout.precision(this->param());
    std::cout << storage.internal_data << std::endl;
  };
};

struct ToString : public SiconosVisitor
{

  std::stringstream sstr;
  std::string answer;

  template<typename T>
  void operator() (const T& storage)
  {
    sstr << storage.internal_data;
    sstr >> answer;
    answer = answer.substr(4, answer.size() - 5); // Remove "[size](" at the beginning of the std::string
    std::string::size_type pos;
    while ((pos = answer.find(",")) != std::string::npos) // Replace "," by " " in the std::string
      answer[pos] = ' ';
  }
};

struct GetValue : public ParamVisitor<unsigned int>
{
  double answer;

  GetValue(unsigned int param) : ParamVisitor(param) {};

  template<typename T>
  void operator() (const T& storage)
  {
    answer = storage.internal_data(this->param());
  };
};

struct GetRValue : public ParamVisitor<unsigned int>
{
  double* answer;

  GetRValue(unsigned int param) : ParamVisitor(param) {};

  template<typename T>
  void operator() (T& storage)
  {
    answer = &(storage.internal_data(this->param()));
  };


};

template<>
void GetRValue::operator()<SparseVectStorage> (SparseVectStorage& storage)
{
  answer = &(storage.internal_data(this->param()).ref());
}

struct SetValue : public ParamVisitor<unsigned int, double>
{

  typedef ParamVisitor<unsigned int, double> base_type;
  SetValue(const base_type::arguments_type& args) :
    base_type(args) {};

  template<typename T>
  void operator() (T& storage) const
  {
    storage.internal_data(this->param1()) = this->param2();
  };
};

struct SetBlock : public ParamVisitor<unsigned int, SiconosVector>
{
  typedef ParamVisitor<unsigned int, SiconosVector> base_type;

  SetBlock(const base_type::arguments_type& args) : base_type(args) {};

  template<typename T>
  void operator() (T& storage)
  {
    unsigned int index = this->param1();
    const SiconosVector& vIn = this->param2();

    if (index >= storage.internal_data.size())
    {
      SiconosVectorException::selfThrow("SiconosVector::setBlock : index ouf of range");
    }

    unsigned int end = vIn.size() + index;
    if (end > storage.internal_data.size())
    {
      SiconosVectorException::selfThrow("SiconosVector::setBlock : invalid ranges");
    }

    noalias(ublas::subrange(storage.internal_data, index, end)) = vIn.dense();

  };
};

/* idem copy, visitor for both parameters */
struct Scal : public ParamVisitor<double>
{

  Scal(const double& value) : ParamVisitor(value) {};

  template<typename T>
  void operator() (T& storage)
  {
    storage.internal_data *= this->param();
  };
};

template<>
void Scal::operator()<DenseVectStorage> (DenseVectStorage& dv)
{
  bindings::scal(this->param(), dv.internal_data);
};


template<typename C>
struct IDot : public ParamVisitor<C>
{

  double answer;

  IDot(const C& p) : ParamVisitor<C>(p) {};

  template<typename T>
  void operator() (const T& storage)
  {
    answer = inner_prod(this->param1().internal_data, storage.internal_data);
  };
};

double inner_prod(DenseVect& v1, DenseVect& v2)
{
  return bindings::dot(v1, v2);
}

struct Dot : public ParamVisitor<SiconosVectorStorage>
{

  double answer;

  Dot(const SiconosVectorStorage& st) : ParamVisitor(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    /* const_cast : VisitorMaker makes visitors on non const references */
    answer = apply_visitor<IDot<T>, double>(const_cast<SiconosVectorStorage&>(this->param()), storage);
  };
};


template<typename C>
struct IOuterProd : public ParamVisitor<C>
{

  SimpleMatrix answer;

  IOuterProd(const C& p) : ParamVisitor<C>(p) {};

  template<typename T>
  void operator() (const T& storage)
  {
    answer = outer_prod(storage.internal_data, this->param1().internal_data);
  };
};


struct OuterProd : public ParamVisitor<SiconosVectorStorage>
{

  SimpleMatrix answer;

  OuterProd(const SiconosVectorStorage& st) : ParamVisitor(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    /* const_cast : VisitorMaker makes visitors on non const references */
    answer = apply_visitor<IOuterProd<T>, SimpleMatrix >(const_cast<SiconosVectorStorage&>(this->param()), storage);
  };
};

template<typename OP, typename S, typename P1=empty, typename P2=empty, typename P3=empty>
struct InternBiOp : public ParamModifier<S, P1, P2, P3>
{
  typedef ParamModifier<S, P1, P2, P3> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
  void operator() (const T& storage)
  {
    static OP op;
    op(this->param1(), this->param2(), this->param3(), storage.internal_data, this->modifiable());
  };

};


template<typename OP, typename S>
struct InternBiOp<OP, S, empty, empty, empty> : public SimpleModifier<S>
{
  typedef SimpleModifier<S> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
    void operator() (const T& storage)
  {
    /* OP modify last argument */
    static OP op;
    op(storage.internal_data, this->modifiable());
  };

};

template<typename OP, typename S, typename P1>
struct InternBiOp<OP, S, P1, empty, empty> : public ParamModifier<S, P1>
{
  typedef ParamModifier<S, P1> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
    void operator() (T& storage)
  {
    static OP op;
    op(this->param1(), storage.internal_data, this->modifiable());
  };

};


template<typename OP, typename S, typename P1, typename P2>
struct InternBiOp<OP, S, P1, P2, empty> : public ParamModifier<S, P1, P2>
{
  typedef ParamModifier<S, P1, P2> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
    void operator() (T& storage)
  {
    static OP op;
    op(this->param1(), this->param2(), storage.internal_data, this->modifiable());
  };

};





/* double pass visitors:
   - store arguments in param
   - visit arg1 and apply internal visitor with visited arg1 as param on param
*/
template<typename OP, typename P1=empty, typename P2=empty, typename P3=empty>
struct BiOperator : public ParamModifier<SiconosVectorStorage, P1, P2, P3>
{

  typedef ParamModifier<SiconosVectorStorage, P1, P2, P3> base_type;
  BiOperator(typename base_type::arguments_type& st) : base_type(st) {};

  template<typename T>
  void operator() (T& storage)
  {
    apply_modifier<InternBiOp<OP, typename T::storage_type, P1, P2, P3> >(this->modifiable(), this->param1(), this->param2(), this->param3(), storage.internal_data);
  };
};

template<typename OP, typename P1, typename P2>
struct BiOperator<OP, P1, P2, empty> :
  public ParamModifier<SiconosVectorStorage, P1, P2>
{

  typedef ParamModifier<SiconosVectorStorage, P1, P2> base_type;
  BiOperator(typename base_type::arguments_type& st) : base_type(st) {};

  template<typename T>
  void operator() (T& storage)
  {
    /* visit and modify modifiable() which is unknown at this level */
    apply_modifier<InternBiOp<OP, typename T::storage_type, P1, P2> >
      (this->modifiable(), this->param1(), this->param2(),
       storage.internal_data);
  };
};

template<typename OP, typename P1>
struct BiOperator<OP, P1, empty, empty> :
  public ParamModifier<SiconosVectorStorage, P1>
{
  typedef ParamModifier<SiconosVectorStorage, P1> base_type;

  BiOperator(typename base_type::arguments_type& st) : base_type(st) {};

  template<typename T>
  void operator() (T& storage)
  {
    /* visit and modify modifiable() which is unknown at this level */
    apply_modifier<InternBiOp<OP, typename T::storage_type, P1> >
      (this->modifiable(), this->param1(), storage.internal_data);
  };
};



template<typename OP>
struct BiOperator<OP, empty, empty, empty> : public ParamVisitor<SiconosVectorStorage>
{
  typedef ParamVisitor<SiconosVectorStorage> base_type;

  /* the storage that is going to be modified is first kept in an argument */
  BiOperator(const typename base_type::arguments_type& st) : base_type(st) {};

  template<typename T>
  void operator() (T& storage)
  {
    /* modify storage.internal_data */
    apply_param_modifier<InternBiOp<OP, typename T::storage_type> >(this->param(), storage.internal_data);
  };
};




struct OpPlus
{
  template<typename T1, typename T2>
  void operator()(const T2& y, T1& x)
  {
    x += y;
  }
};

struct OpMinus
{
  template<typename T1, typename T2>
  void operator()(const T2& y, T1& x)
  {
    x -= y;
  }
};

struct OpCopy
{
  template<typename T1, typename T2>
  void operator()(const T2& y, T1& x)
  {
    x.resize(y.size());
    x = y;
  };
};

template<>
void OpCopy::operator()(const DenseVect& dv2, DenseVect& dv1)
{
  dv1.resize(dv2.size());
  bindings::copy(dv2, dv1);
}

struct OpAxpy
{
  template<typename T1, typename T2>
  void operator()(const double& param, const T2& x, T1& y)
  {
    T2 tmp = param*x;
    y += tmp;
  }
};

struct OpAxpby
{
  template<typename T1, typename T2>
  void operator()(const double& a, const double& b, const T2& x, T1& y)
  {
    y *= b;
    y += a * x;
  }
};

template<>
void OpAxpby::operator()<DenseVect, DenseVect>(const double& a,
                                               const double& b,
                                               const DenseVect& x,
                                               DenseVect& y)
{
  bindings::scal(b, y);
  bindings::axpy(a, x, y);
}

template<>
void OpAxpby::operator()<SparseVect, SparseVect>(const double& a,
                                                 const double& b,
                                                 const SparseVect& x,
                                                 SparseVect& y)
{
  y *= b;
  if (&y != &x)
  {
    noalias(y) += a*x;
  }
  else
  {
    y += a*x;
  }
}


struct OpAddBlock
{
  template<typename T1, typename T2>
  void operator() (unsigned int param1, const T1& storage1, T2& storage2)
  {
    unsigned int index = param1;
    unsigned int end = storage2.size();
    if ((index + end) > storage1.size())
    {
      SiconosVectorException::selfThrow("SiconosVector::addBlock : invalid ranges");
    }

    noalias(ublas::subrange(storage2, index, index + end)) += storage1;
  };
};

struct OpSubBlock
{
  template<typename T1, typename T2>
  void operator() (unsigned int param1, const T1& storage1, T2& storage2)
  {
    unsigned int index = param1;
    unsigned int end = storage2.size();
    if ((index + end) > storage1.size())
    {
      SiconosVectorException::selfThrow("SiconosVector::subBlock : invalid ranges");
    }

    noalias(ublas::subrange(storage2, index, index + end)) -= storage1;
  };
};

struct OpToBlock
{
  template <typename T1, typename T2>
  void operator() (unsigned int param1, unsigned int param2, unsigned int param3, const T1& storage1, T2& storage2)
  {
    T2& vOut = storage2;
    unsigned int sizeB = param1;
    unsigned int startIn = param2;
    unsigned int startOut = param3;

    // To copy a subBlock of the vector (from position startIn to startIn+sizeB) into vOut (from pos. startOut to startOut+sizeB).
    // Check dim ...
    assert(startIn < storage1.size() && "vector toBlock(v1,v2,...): start position in input vector is out of range.");

    assert(startOut < vOut.size() && "vector toBlock(v1,v2,...): start position in output vector is out of range.");

    assert(startIn + sizeB <= storage1.size() && "vector toBlock(v1,v2,...): end position in input vector is out of range.");
    assert(startOut + sizeB <= vOut.size() && "vector toBlock(v1,v2,...): end position in output vector is out of range.");

    unsigned int endOut = startOut + sizeB;

    noalias(ublas::subrange(vOut, startOut, endOut)) = ublas::subrange(storage1, startIn, startIn + sizeB);
  };
};


struct OpSubscal
{
  template<typename T1, typename T2>
  void operator() (const double& param1, const Index& param2, bool param3, const T1& storage1, T2& storage2)
  {
    const double& a = param1;
    const Index& coord = param2;
    bool init = param3;

    unsigned int dimX = coord[1] - coord[0];
    unsigned int dimY = coord[3] - coord[2];

    const T1& x = storage1;
    T2& y = storage2;

    if (dimY != dimX)
    {
      SiconosVectorException::selfThrow("subscal(a,x,y,...) error: inconsistent sizes between (sub)x and (sub)y.");
    }
    if (dimY > y.size() || dimX > x.size())
    {
      SiconosVectorException::selfThrow("subscal(a,x,y,...) error: input index too large.");
    }
    if((void*)&x == (void*)&y)
    {
      ublas::vector_range<T2> subY(y, ublas::range(coord[2], coord[3]));
      if (coord[0] == coord[2])
      {
        if (init)
          subY *= a;
        else
          subY *= (1.0 + a);
      }
      else
      {
        const ublas::vector_range<const T1> subX(x, ublas::range(coord[0], coord[1]));
        if (init)
          subY = a * subX;
        else
          subY += a * subX;
      }
    }
    else
    {
      const ublas::vector_range<const T1> subX(x, ublas::range(coord[0], coord[1]));
      ublas::vector_range<T2> subY(y, ublas::range(coord[2], coord[3]));
      if (init)
      {
        noalias(subY) = a * subX;
      }
      else
      {
        noalias(subY) += a * subX;
      }
    }
  }
};


struct Plus : public BiOperator<OpPlus>
{
  typedef BiOperator<OpPlus> base_type;
  Plus(const base_type::arguments_type& args) : base_type(args) {};
};
struct Minus : public BiOperator<OpMinus>
{
  typedef BiOperator<OpMinus> base_type;
  Minus(const base_type::arguments_type& args) : base_type(args) {};
};
struct Copy : public BiOperator<OpCopy>
{
  typedef BiOperator<OpCopy> base_type;
  Copy(const base_type::arguments_type& args) : base_type(args) {};
};
struct AddBlock : public BiOperator<OpAddBlock, unsigned int>
{
  typedef BiOperator<OpAddBlock, unsigned int> base_type;
  AddBlock(base_type::arguments_type& args) : base_type(args) {};
};
struct SubBlock : public BiOperator<OpSubBlock, unsigned int>
{
  typedef BiOperator<OpSubBlock, unsigned int> base_type;
  SubBlock(base_type::arguments_type& args) : base_type(args) {};
};
struct Axpy : public BiOperator<OpAxpy, double>
{
  typedef BiOperator<OpAxpy, double> base_type;
  Axpy(base_type::arguments_type& args) : base_type(args) {};
};

struct Axpby : public BiOperator<OpAxpby, double, double>
{
  typedef BiOperator<OpAxpby, double, double> base_type;
  Axpby(base_type::arguments_type& args) : base_type(args) {};
};

struct Subscal : public BiOperator<OpSubscal, double, const Index, bool>
{
  typedef BiOperator<OpSubscal, double, const Index, bool> base_type;
  Subscal(base_type::arguments_type& args) : base_type(args) {};
};
struct ToBlock : public BiOperator<OpToBlock, unsigned int, unsigned int, unsigned int>
{
  typedef BiOperator<OpToBlock, unsigned int, unsigned int, unsigned int> base_type;
  ToBlock(base_type::arguments_type& args) : base_type(args) {};
};

struct StorageAllocator : public ParamVisitor<unsigned int>
{
  SiconosVectorStorage* answer;

  StorageAllocator(unsigned int size) : ParamVisitor(size) {};

  template<typename T>
  void operator ()(const T&)
  {
    answer = new T(this->param());
  }
};

// =================================================
//                CONSTRUCTORS
// =================================================

// Default
SiconosVector::SiconosVector() : _storage(new DenseVectStorage()) {};

// parameters: dimension and type.
SiconosVector::SiconosVector(unsigned row, Siconos::UBLAS_TYPE type)
{
  if (type == Siconos::SPARSE)
  {
    _storage = new SparseVectStorage(row);
  }
  else if (type == Siconos::DENSE)
  {
    _storage = new DenseVectStorage(row);
  }
  else
  {
    SiconosVectorException::selfThrow("SiconosVector::constructor(Siconos::UBLAS_TYPE, unsigned int) failed, invalid type given");
  }
}

// parameters: dimension, default value for all components and type.
SiconosVector::SiconosVector(unsigned row, double val, Siconos::UBLAS_TYPE type)
{
  if (type == Siconos::SPARSE)
  {
    _storage = new SparseVectStorage(row);
    fill(val);
  }
  else if (type == Siconos::DENSE)
  {
    _storage = new DenseVectStorage(row, val);
  }
  else
  {
    SiconosVectorException::selfThrow("SiconosVector::constructor(Siconos::UBLAS_TYPE, unsigned int) : invalid type given");
  }
}

// parameters: a vector (stl) of double and the type.
SiconosVector::SiconosVector(const std::vector<double>& v, Siconos::UBLAS_TYPE typ)
{
  if (typ != Siconos::DENSE)
    SiconosVectorException::selfThrow("SiconosVector::constructor(Siconos::UBLAS_TYPE, std::vector<double>, unsigned int) : invalid type given");

 _storage = new DenseVectStorage(v.size());
 std::copy(v.begin(), v.end(), static_cast<DenseVectStorage*>(_storage)->internal_data.begin());
}

// Copy
SiconosVector::SiconosVector(const SiconosVector &svect) :
  std11::enable_shared_from_this<SiconosVector>()
{
  this->_storage = apply_visitor<StorageAllocator, SiconosVectorStorage*>(storage(svect), svect.size());
  apply_modifier<Copy>(storage(*this), storage(svect));
}

// Copy from BlockVector
SiconosVector::SiconosVector(const BlockVector & vIn) : std11::enable_shared_from_this<SiconosVector>()
{
  this->_storage = apply_visitor<StorageAllocator, SiconosVectorStorage*>(storage(**vIn.begin()), vIn.size());
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    setBlock(pos, **it);
    pos += (*it)->size();
  }

}

SiconosVector::SiconosVector(const DenseVect& m)
{
  this->_storage = new DenseVectStorage(m.size());
  noalias(this->dense()) = m;

}

SiconosVector::SiconosVector(const SparseVect& m)
{
  this->_storage = new SparseVectStorage(m.size());
  noalias(this->sparse()) = m;
}

SiconosVector::SiconosVector(const std::string &file, bool ascii)
{
  this->_storage = new DenseVectStorage();
  if (ascii)
  {
    ioVector::read(file, *this, ioVector::ASCII_IN);
   }
  else
  {
    ioVector::read(file, *this, ioVector::BINARY_IN);
  }
}

SiconosVector::SiconosVector(const SiconosVector& v1, const SiconosVector& v2)
{
  unsigned int size1 = v1.size();

  if (Type::value(storage(v1)) == Type::value(storage(v2)))
  {
    this->_storage = apply_visitor<StorageAllocator, SiconosVectorStorage*>(storage(v1), v1.size() + v2.size()) ;
  }
  else
  {
    SiconosVectorException::selfThrow("SiconosVector::SiconosVector :: mixed dense and sparse vector detected");
  }
  setBlock(0, v1);
  setBlock(size1, v2);
}

SiconosVector::~SiconosVector()
{
  if (_storage) delete(_storage);
}


// =================================================
//        get Ublas component (dense or sparse)
// =================================================

unsigned int SiconosVector::num() const
{
  if (Type::value(*_storage) == Type::DenseVectStorage)
  {
    return 1;
  }
  else
    if (Type::value(*_storage) == Type::SparseVectStorage)
    {
      return 4;
    }
    else
    {
      return 0;
    }
}

DenseVect& SiconosVector::dense(unsigned int) const
{
  if (Type::value(*_storage) != Type::DenseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::Densevect(unsigned int) : cannot get dense storage.");
  }
  return static_cast<DenseVectStorage*>(_storage)->internal_data;
}

SparseVect& SiconosVector::sparse(unsigned int)const
{
  if (Type::value(*_storage) != Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::Sparsevect(unsigned int) : cannot get sparse storage.");
  }
  return static_cast<SparseVectStorage*>(_storage)->internal_data;
}

double* SiconosVector::getArray() const
{
  if (Type::value(*_storage) != Type::DenseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::getArray() : cannot get array for this kind of vector.");
  }
  return &(((static_cast<DenseVectStorage*>(_storage)->internal_data).data())[0]);
}

// ===========================
//       fill vector
// ===========================

void SiconosVector::zero()
{
  apply_modifier<Zero>(storage(*this));
}

void SiconosVector::setVector(unsigned int , const SiconosVector& newV)
{
  if (newV.size() != size())
    SiconosVectorException::selfThrow("SiconosVector::setVector(num,v), unconsistent sizes.");

  *this = newV ;
}

void SiconosVector::fill(const double& value)
{
  apply_modifier<Fill>(storage(*this), value);
}

//=======================
// set vector dimension
//=======================

void SiconosVector::resize(unsigned int n, bool preserve)
{
  apply_modifier<Resize>(storage(*this), n);
}

//=======================
//       get norm
//=======================

double SiconosVector::normInf() const
{
  return apply_visitor<NormInf, double>(storage(*this));
}

double SiconosVector::norm2() const
{
  return apply_visitor<Norm2, double>(storage(*this));
}
//======================================
// get sum of all elements of the vector
//=====================================
double SiconosVector::vector_sum() const
{
  return apply_visitor<Sum, double>(storage(*this));
}

//=====================
// screen display
//=====================

void SiconosVector::display(unsigned int n)const
{
  apply_visitor<Display>(storage(*this), n);
}

//============================
// Convert vector to a std::string
//============================

const std::string SiconosVector::toString() const
{
  return apply_visitor<ToString, std::string>(storage(*this));
}

//=============================
// Elements access (get or set)
//=============================

double SiconosVector::getValue(unsigned int row) const
{
  return apply_visitor<GetValue, unsigned int>(storage(*this), row);
}

void SiconosVector::setValue(unsigned int row, const double value)
{
  return apply_modifier<SetValue>(storage(*this), row, value);
}

double& SiconosVector::operator()(unsigned int row)
{
  return *apply_modifier<GetRValue, double* >(storage(*this), row);
};

double SiconosVector::operator()(unsigned int row) const
{
  return getValue(row);
}

// //============================================
// // Access (get or set) to blocks of elements
// //============================================

void SiconosVector::setBlock(unsigned int index, const SiconosVector& vIn)
{
  apply_modifier<SetBlock>(storage(*this), index, vIn);
}

void SiconosVector::toBlock(SiconosVector& vOut, unsigned int sizeB, unsigned int startIn, unsigned int startOut) const
{
  apply_modifier<ToBlock>(storage(*this), sizeB, startIn, startOut, storage(vOut));
}

void SiconosVector::addBlock(unsigned int index, const SiconosVector& vIn)
{
  apply_modifier<AddBlock>(storage(vIn), index, storage(*this));
}

void SiconosVector::subBlock(unsigned int index, const SiconosVector& vIn)
{
  apply_modifier<SubBlock>(storage(vIn), index, storage(*this));
}

// //===============
// //  Assignment
// //===============

SiconosVector& SiconosVector::operator = (const SiconosVector& vIn)
{
  if (&vIn == this)
  {
    return *this; // auto-assignment.
  }
  else
  {
    apply_modifier<Copy>(storage(*this), storage(vIn));
    return *this;
  }
}

SiconosVector& SiconosVector::operator = (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    setBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}


SiconosVector& SiconosVector::operator = (const DenseVect& d)
{
  if (Type::value(storage(*this)) != Type::DenseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : current vector is not dense.");
  }
  if (this->size() != d.size())
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = DenseVect : inconsistent size.");
  }
  bindings::copy(d, this->dense());
  return *this;
}

SiconosVector& SiconosVector::operator = (const SparseVect& sp)
{
  if (Type::value(storage(*this)) != Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : current vector is not sparse.");
  }
  if (this->size() != sp.size())
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = SparseVect : inconsistent size.");
  }

  noalias(this->sparse()) = sp;

  return *this;
}

SiconosVector& SiconosVector::operator = (const double* d)
{
  if (Type::value(storage(*this)) == Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::operator = double* : forbidden: the current vector is not dense.");
  }
  bindings::detail::copy(this->size(), d, 1, getArray(), 1);
  return *this;
}

unsigned SiconosVector::copyData(double* data) const
{
  if (Type::value(storage(*this)) == Type::SparseVectStorage)
  {
    SiconosVectorException::selfThrow("SiconosVector::copyData : forbidden: the current vector is not dense.");
  }
  unsigned size = this->size();
  bindings::detail::copy(size, getArray(), 1, data, 1);
  return size;
}


//=================================
// Op. and assignment (+=, -= ... )
//=================================

SiconosVector& SiconosVector::operator += (const SiconosVector& vIn)
{
  apply_modifier<Plus>(storage(*this), storage(vIn));
  return *this;
}
SiconosVector& SiconosVector::operator += (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    addBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}

SiconosVector& SiconosVector::operator -= (const SiconosVector& vIn)
{
  apply_modifier<Minus>(storage(*this), storage(vIn));
  return *this;
}

SiconosVector& SiconosVector::operator -= (const BlockVector& vIn)
{
  VectorOfVectors::const_iterator it;
  unsigned int pos = 0;
  for (it = vIn.begin(); it != vIn.end(); ++it)
  {
    subBlock(pos, **it);
    pos += (*it)->size();
  }
  return *this;
}


//===============
// Comparison
//===============

bool operator == (const SiconosVector &m, const SiconosVector &x)
{
  DEBUG_PRINTF("norm = %12.8e \n", (m - x).normInf() );
  DEBUG_PRINTF("std::numeric_limits<double>::epsilon() = %12.8e \n", std::numeric_limits<double>::epsilon() );
  DEBUG_EXPR(std::cout << std::boolalpha << ( (m - x).normInf() <= std::numeric_limits<double>::epsilon()) <<std::endl;);
  return ((m - x).normInf() <= std::numeric_limits<double>::epsilon());
}

//==================
// y = scalar * x
//==================

SiconosVector operator * (const  SiconosVector&m, double d)
{
  SiconosVector tmp = m;
  apply_modifier<Scal>(storage(tmp), d);
  return tmp;
}

SiconosVector operator * (double d, const  SiconosVector&m)
{
  SiconosVector tmp = m;
  apply_modifier<Scal>(storage(tmp), d);
  return tmp;
}

SiconosVector operator / (const SiconosVector &m, double d)
{
  return m * (1.0/d);
}

//====================
//  Vectors addition
//====================

SiconosVector operator + (const  SiconosVector& x, const  SiconosVector& y)
{
  SiconosVector tmp = x;
  tmp += y;
  return tmp;
}

void add(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  apply_modifier<Copy>(storage(z), storage(x));
  z += y;
}

//======================
//  Vectors subtraction
//======================

SiconosVector operator - (const  SiconosVector& x, const  SiconosVector& y)
{
  SiconosVector tmp = x;
  tmp -= y;
  return tmp;
}

void sub(const SiconosVector& x, const SiconosVector& y, SiconosVector& z)
{
  apply_modifier<Copy>(storage(z), storage(x));
  z -= y;
}




void axpby(double a, const SiconosVector& x, double b, SiconosVector& y)
{
  apply_modifier<Axpby>(storage(x), a, b, storage(y));
}

void axpy(double a, const SiconosVector& x, SiconosVector& y)
{
  apply_modifier<Axpy>(storage(x), a, storage(y));
}


double inner_prod(const SiconosVector &x, const SiconosVector &m)
{
  return apply_visitor<Dot, double>(storage(x), storage(m));
}

//// outer_prod(v,w) = trans(v)*w
SimpleMatrix outer_prod(const SiconosVector &x, const SiconosVector& m)
{
  return apply_visitor<OuterProd, SimpleMatrix >(storage(x), storage(m));
}

void scal(double a, const SiconosVector & x, SiconosVector & y, bool init)
{
  if(init)
  {
    apply_modifier<Copy>(storage(y), storage(x));
    apply_modifier<Scal>(storage(y), a);
  }
  else
  {
    apply_modifier<Axpy>(storage(x), a, storage(y));
  }
}

void subscal(double a, const SiconosVector & x, SiconosVector & y, const Index& coord, bool init)
{
  apply_modifier<Subscal>(storage(x), a, coord, init, storage(y));
}
void cross_product(const SiconosVector& V1, const SiconosVector& V2, SiconosVector& VOUT)
{
  if (V1.size() != 3 || V2.size() != 3 || VOUT.size() != 3)
    SiconosVectorException::selfThrow("SiconosVector::cross_product allowed only with dim 3.");

  double aux = V1.getValue(1) * V2.getValue(2) - V1.getValue(2) * V2.getValue(1);
  VOUT.setValue(0, aux);

  aux = V1.getValue(2) * V2.getValue(0) - V1.getValue(0) * V2.getValue(2);
  VOUT.setValue(1, aux);

  aux = V1.getValue(0) * V2.getValue(1) - V1.getValue(1) * V2.getValue(0);
  VOUT.setValue(2, aux);

}

//

void abs_wise(const SiconosVector& V, SiconosVector& Vabs)
{
  for (unsigned int it = 0; it < V.size(); ++it)
  {
    Vabs.setValue(it, std::abs(V.getValue(it)));
  };
}

//

void getMax(const SiconosVector& V, double& maxvalue, unsigned int& idmax)
{
  maxvalue = V.getValue(0);
  idmax = 0;
  for (unsigned int it = 1; it < V.size(); ++it)
  {
    if (V.getValue(it) > maxvalue)
    {
      maxvalue = V.getValue(it);
      idmax = it;
    };
  };
}

//

void getMin(const SiconosVector& V, double& minvalue, unsigned int& idmin)
{
  minvalue = V.getValue(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it)
  {
    if (V.getValue(it) < minvalue)
    {
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}


struct exp_op { double operator() (double d) const { return std::exp(d); } };

void SiconosVector::exp_in_place()
{
  // struct exp_op { double operator() (double d) const { return std::exp(d); } };
  // assert(num() == 1);
  // std::transform(vect.Dense->begin(), vect.Dense->end(), vect.Dense->begin(), exp_op);
}

void SiconosVector::exp(SiconosVector& input)
{
  // assert(num() == 1 && input.num()==1);
  // std::transform(input.dense()->begin(), input.dense()->end(), vect.Dense->begin(), exp_op);
}


//
/*
SiconosVector abs_wise(const SiconosVector& V){
  SiconosVector Vabs(V.size());
  for (int it = 0; it < V.size(); ++it){
    Vabs.setValue(it,std::abs(V.getValue(it)));
  };
  return Vabs;
}
//
void getMin(const SiconosVector& V, double& minvalue, unsigned int& idmin){
  minvalue = V.getValue(0);
  idmin = 0;
  for (unsigned int it = 1; it < V.size(); ++it){
    if (V.getValue(it) < minvalue){
      minvalue = V.getValue(it);
      idmin = it;
    };
  };
}
*/
void setBlock(const SiconosVector& vIn, SP::SiconosVector vOut, unsigned int sizeB,
              unsigned int startIn, unsigned int startOut)
{
  apply_modifier<ToBlock>(const_cast<SiconosVectorStorage&>(storage(vIn)), sizeB, startIn, startOut, const_cast<SiconosVectorStorage&>(storage(*vOut)));
}

unsigned int SiconosVector::size(void) const
{
  return apply_visitor<Size, unsigned int>(storage(*this));
}

SiconosVector& operator *= (SiconosVector& v, const double& s)
{

  apply_modifier<Scal>(storage(v), s);
  return v;
}


SiconosVector& operator /= (SiconosVector& v, const double& s)
{
  apply_modifier<Scal>(storage(v), 1.0/s);
  return v;
}
