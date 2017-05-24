#ifndef SiconosVectorStorage_hpp
#define SiconosVectorStorage_hpp

#include "SiconosAlgebraTypeDef.hpp"
#include "SiconosVisitor.hpp"
#include "SiconosVectorException.hpp"

#include <boost/array.hpp>
#include <vector>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>

#include <boost/move/unique_ptr.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/utility/enable_if.hpp>

class SiconosVectorStorage
{
public:
  ACCEPT_STD_VISITORS();

  virtual ~SiconosVectorStorage() {};
  SiconosVectorStorage() {};
private:
  SiconosVectorStorage(SiconosVectorStorage&) {};

};

struct DenseVectStorage : public SiconosVectorStorage
{
  typedef DenseVect storage_type;
  DenseVect internal_data;
  DenseVectStorage();
  DenseVectStorage(size_t r);

  DenseVectStorage(size_t r, double val);

  virtual ~DenseVectStorage() {};

  ACCEPT_STD_VISITORS();
};

template<size_t N>
struct BoundedVectStorage : public SiconosVectorStorage
{
  typedef ublas::vector<double, ublas::bounded_array<double, N> > storage_type;
  storage_type internal_data;
  BoundedVectStorage() : internal_data() {};
  BoundedVectStorage(size_t r) : internal_data(r) {};

  BoundedVectStorage(size_t r, double val) : internal_data(r, val) {};

  virtual ~BoundedVectStorage() {};

  ACCEPT_STD_VISITORS();
};

/* issue with forwarding typedefs */
#define BOUNDED_VECT_STORAGE_INSTANCE(N) \
  struct BoundedVectStorage##N : public BoundedVectStorage<N>                       \
  {                                                                                 \
    BoundedVectStorage##N() : BoundedVectStorage<N>() {};                           \
    BoundedVectStorage##N(size_t r) : BoundedVectStorage<N>(r) {};                  \
    BoundedVectStorage##N(size_t r, double val) : BoundedVectStorage<N>(r, val) {}; \
                                                                                    \
   ACCEPT_STD_VISITORS();                                                           \
  }

BOUNDED_VECT_STORAGE_INSTANCE(3);
BOUNDED_VECT_STORAGE_INSTANCE(4);
BOUNDED_VECT_STORAGE_INSTANCE(6);
BOUNDED_VECT_STORAGE_INSTANCE(7);

struct SparseVectStorage : public SiconosVectorStorage
{
  typedef SparseVect storage_type;
  SparseVect internal_data;
  SparseVectStorage();
  SparseVectStorage(size_t r);

  virtual ~SparseVectStorage() {};

  ACCEPT_STD_VISITORS();
};

#ifndef SWIG

#include <boost/type_traits.hpp>
#include <boost/mpl/if.hpp>

/* to be forwarded as a friend class */
struct Storage
{
  template<typename T, typename V>
  struct forward_const
  {
    typedef typename boost::mpl::if_<boost::is_const<T>, const V, V>::type type;
  };

  template<typename T>
  static typename forward_const<T, SiconosVectorStorage>::type& get(T& v)
  {
    return *v._storage;
  }
};


template<typename T>
static typename Storage::forward_const<T, SiconosVectorStorage>::type& storage(T& v)
{
  return Storage::get(v);
};

#undef VISITOR_CLASSES
#define VISITOR_CLASSES()                       \
  REGISTER(DenseVectStorage)                    \
  REGISTER(BoundedVectStorage3)                 \
  REGISTER(BoundedVectStorage4)                 \
  REGISTER(BoundedVectStorage6)                 \
  REGISTER(BoundedVectStorage7)                 \
  REGISTER(SparseVectStorage)

#include <VisitorMaker.hpp>

#include <iostream>
#include <typeinfo>

#include <cmath>        // std::exp(double)
#include <algorithm>    // std::transform


using namespace Experimental;

template<typename T>
struct const_argument
{
  typedef typename boost::mpl::if_<boost::is_integral<T>, const T, const T&>::type type;
};

template<typename T>
struct ref_argument
{
  typedef typename boost::mpl::if_<boost::is_integral<T>, T, T&>::type type;
};

template<typename T>
struct build_argument
{
  typedef typename boost::remove_const<typename boost::remove_reference<T>::type>::type type;
};

template<typename T, typename R>
struct enable_if_is_number
{
  typedef typename boost::enable_if<boost::mpl::or_<boost::is_integral<T>, boost::is_floating_point<T> >, R>::type type;
};

template<typename T, typename R>
struct disable_if_is_number
{
  typedef typename boost::disable_if<boost::mpl::or_<boost::is_integral<T>, boost::is_floating_point<T>, boost::is_pointer<T> >, R>::type type;
};

/* a visitor with some const arguments */
template<typename P1, typename P2=empty, typename P3=empty, typename P4=empty>
struct ParamVisitor : public SiconosVisitor
{
  typedef boost::tuple<P1, P2, P3, P4> arguments_type;

  const arguments_type& _param;

  ParamVisitor(const arguments_type& args) : _param(args) {};

  virtual ~ParamVisitor() {};

  P1 param1() const { return boost::get<0>(_param); };
  P2 param2() const { return boost::get<1>(_param); };
  P3 param3() const { return boost::get<2>(_param); };
  P4 param4() const { return boost::get<3>(_param); };
};


template<>
struct ParamVisitor<empty, empty, empty, empty> : public SiconosVisitor
{
  virtual ~ParamVisitor() {};
};

template<typename P>
struct ParamVisitor<P, empty, empty, empty> : public SiconosVisitor
{
  typedef P arguments_type;
  typename ref_argument<arguments_type>::type _param;

  ParamVisitor(typename ref_argument<arguments_type>::type a) : _param(a) {};

  virtual ~ParamVisitor() {};

  typename ref_argument<arguments_type>::type param() { return _param; };
};

template<typename P1, typename P2>
struct ParamVisitor<P1, P2, empty, empty> : public SiconosVisitor
{
  typedef boost::tuple<P1, P2> arguments_type;
  const arguments_type& _param;

  ParamVisitor(const arguments_type& args) : _param(args) {};

  virtual ~ParamVisitor() {};

  P1 param1() const { return boost::get<0>(_param); };
  P2 param2() const { return boost::get<1>(_param); };

};

template<typename P1, typename P2, typename P3>
struct ParamVisitor<P1, P2, P3, empty> : public SiconosVisitor
{
  typedef boost::tuple<P1, P2, P3> arguments_type;
  const arguments_type& _param;

  ParamVisitor(const arguments_type& args) : _param(args) {};

  virtual ~ParamVisitor() {};

  P1 param1() const { return boost::get<0>(_param); };
  P2 param2() const { return boost::get<1>(_param); };
  P3 param3() const { return boost::get<2>(_param); };

};

#undef REGISTER
#define REGISTER(X) X,
template <typename V, typename IsConst>
struct make_visitor : public Visitor<Classes< VISITOR_CLASSES() empty>, V, typename boost::mpl::not_<IsConst> >::Make
{

  typedef typename Visitor<Classes< VISITOR_CLASSES() empty>, V, typename boost::mpl::not_<IsConst> >::Make base_type;

  make_visitor()  {};

  template<typename T>
  make_visitor(T& args) :
    base_type(args) {};

  make_visitor(unsigned int args) :
    base_type(args) {};

};

template<typename V, typename T>
void apply_visitor(V& visitor, T& obj)
{
  obj.accept(visitor);
};

template<typename V, typename T>
void apply_visitor(T& obj)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor = make_visitor<V, typename boost::is_const<T>::type>();
  obj.accept(visitor);
};

template<typename V, typename T, typename P>
void apply_visitor(T& obj, P& param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
};

template<typename V, typename T, typename P>
void apply_visitor(T& obj, const P& param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
};

template<typename V, typename T>
void apply_visitor(T& obj, unsigned int param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
};

template<typename V, typename T>
void apply_visitor(T& obj, int param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
};

template<typename V, typename T>
void apply_visitor(T& obj, double param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1, typename P2>
void apply_visitor(T& obj, P1 param1, P2& param2)
{
  typedef typename V::arguments_type P;
  P args(param1, param2);
  make_visitor<V, typename boost::is_const<T>::type> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1>
void apply_visitor(T& obj, P1 param1, double param2)
{
  typedef typename V::arguments_type P;
  P args(param1, param2);
  make_visitor<V, typename boost::is_const<T>::type> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1, typename P2, typename P3>
void apply_visitor(T& obj, P1 param1, P2 param2, P3& param3)
{
  typedef typename V::arguments_type P;
  P args(param1, param2, param3);
  make_visitor<V, typename boost::is_const<T>::type> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1, typename P2, typename P3, typename P4>
void apply_visitor(T& obj, P1& param1, P2& param2, P3& param3, P4& param4)
{
  typedef typename V::arguments_type P;
  P args(param1, param2, param3, param4);
  make_visitor<V, typename boost::is_const<T>::type> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename T, typename P1, typename P2, typename P3, typename P4>
void apply_visitor(T& obj, P1& param1, P2& param2, const P3& param3, P4& param4)
{
  typedef typename V::arguments_type P;
  P args(param1, param2, param3, param4);
  make_visitor<V, typename boost::is_const<T>::type> visitor(args);
  obj.accept(visitor);
};

template<typename V, typename R, typename T>
R apply_visitor(T& obj)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor;
  obj.accept(visitor);
  return visitor.answer;
};

template<typename V, typename R, typename T>
R apply_visitor(T& obj, unsigned int param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
  return visitor.answer;
};

template<typename V, typename R, typename T>
R apply_visitor(T& obj, const int param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
  return visitor.answer;
};

template<typename V, typename R, typename T>
R apply_visitor(T& obj, double param)
{
  make_visitor<V, typename boost::is_const<T>::type> visitor(param);
  obj.accept(visitor);
  return visitor.answer;
};

template<typename V, typename R, typename P, typename T>
R apply_visitor(T& obj, P& param)
{
  make_visitor<V, typename boost::is_const<T>::type > visitor(param);
  obj.accept(visitor);
  return visitor.answer;
};

#endif
/* --- */
#endif
