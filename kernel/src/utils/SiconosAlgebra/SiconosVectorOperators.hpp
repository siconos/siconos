#ifndef SiconosVectorOperators_hpp
#define SiconosVectorOperators_hpp
#include "SiconosVectorStorage.hpp"

namespace bindings = boost::numeric::bindings::blas;

struct Resize : public ParamVisitor<const size_t>
{

  Resize(size_t size) : ParamVisitor(size) {};

  template<typename T>
  void operator() (T& v)
  {
    assert (this->param() >= 0);
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

struct Sum : public SiconosVisitor
{
  double answer;

  template<typename T>
  void operator () (const T& storage)
  {
    answer = sum(storage.internal_data);
  };
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

struct GetValue : public ParamVisitor<const unsigned int>
{
  double answer;

  GetValue(const unsigned int param) : ParamVisitor(param) {};

  template<typename T>
  void operator() (const T& storage)
  {
    assert (this->param() >= 0);
    assert (this->param() < storage.internal_data.size());

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
    assert (this->param() >= 0);
    assert (this->param() <= storage.internal_data.size());

    answer = &(storage.internal_data(this->param()));
  };


};

template<>
void GetRValue::operator()<SparseVectStorage> (SparseVectStorage& storage)
{
  assert (this->param() <= storage.internal_data.size());

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
    assert (this->param1() >= 0);
    assert (this->param1() <= storage.internal_data.size());

    storage.internal_data(this->param1()) = this->param2();
  };
};

struct OpSetBlock
{

  template<typename T1, typename T2>
  void operator() (const T1& storage1, unsigned int param1, T2& storage2)
  {
    unsigned int index = param1;

    if (index >= storage2.size())
    {
      SiconosVectorException::selfThrow("SiconosVector::setBlock : index ouf of range");
    }

    unsigned int end = storage1.size() + index;
    if (end > storage2.size())
    {
      SiconosVectorException::selfThrow("SiconosVector::setBlock : invalid ranges");
    }

    noalias(ublas::subrange(storage2, index, end)) = storage1;

  };
};

/* idem copy, visitor for both parameters */
struct Scal : public ParamVisitor<const double>
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
struct IDot : public ParamVisitor<const C&>
{

  double answer;

  IDot(const C& p) : ParamVisitor<const C&>(p) {};

  template<typename T>
  void operator() (const T& storage)
  {
    answer = inner_prod(this->param().internal_data, storage.internal_data);
  };
};

double inner_prod(DenseVect& v1, DenseVect& v2)
{
  return bindings::dot(v1, v2);
}

struct Dot : public ParamVisitor<const SiconosVectorStorage&>
{

  double answer;

  Dot(const SiconosVectorStorage& st) : ParamVisitor(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    answer = apply_visitor<IDot<T>, double>(this->param(), storage);
  };
};

template<typename OP, typename P1=empty, typename P2=empty, typename P3=empty, typename P4=empty>
struct InternBiOp : public ParamVisitor<P1, P2, P3, P4>
{
  typedef ParamVisitor<P1, P2, P3, P4> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
  void operator() (T& storage)
  {
    static OP op;
    op(this->param4(), this->param1(), this->param2(), this->param3(), storage.internal_data);
  };

};


template<typename OP, typename S>
struct InternBiOp<OP, S, empty, empty, empty> : public ParamVisitor<S>
{
  typedef ParamVisitor<S> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
    void operator() (T& storage)
  {
    /* OP modify last argument */
    static OP op;
    op(this->param(), storage.internal_data);
  };

};

template<typename OP, typename P1, typename P2>
struct InternBiOp<OP, P1, P2, empty, empty> : public ParamVisitor<P1, P2>
{
  typedef ParamVisitor<P1, P2> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
    void operator() (T& storage)
  {
    static OP op;
    op(this->param2(), this->param1(), storage.internal_data);
  };

};


template<typename OP, typename P1, typename P2, typename P3>
struct InternBiOp<OP, P1, P2, P3, empty> : public ParamVisitor<P1, P2, P3>
{
  typedef ParamVisitor<P1, P2, P3> base_type;

  InternBiOp(typename base_type::arguments_type& p) : base_type(p) {};

  template<typename T>
    void operator() (T& storage)
  {
    static OP op;
    op(this->param3(), this->param1(), this->param2(), storage.internal_data);
  };

};





/* double pass visitors:
   - store arguments in param
   - visit arg1 and apply internal visitor with visited arg1 as param on param
*/
template<typename OP, typename P1=empty, typename P2=empty, typename P3=empty>
struct BiOperator : public ParamVisitor<P1, P2, P3, SiconosVectorStorage&>
{

  typedef ParamVisitor<P1, P2, P3, SiconosVectorStorage&> base_type;
  BiOperator(const typename base_type::arguments_type& st) : base_type(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    apply_visitor<InternBiOp<OP, P1, P2, P3, const typename T::storage_type&> >(this->param4(), this->param1(), this->param2(), this->param3(), storage.internal_data);
  };
};

template<typename OP, typename P1, typename P2>
struct BiOperator<OP, P1, P2, empty> :
  public ParamVisitor<P1, P2, SiconosVectorStorage&>
{

  typedef ParamVisitor<P1, P2, SiconosVectorStorage&> base_type;
  BiOperator(const typename base_type::arguments_type& st) : base_type(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    apply_visitor<InternBiOp<OP, P1, P2, const typename T::storage_type&> >
      (this->param3(), this->param1(), this->param2(), storage.internal_data);
  };
};

template<typename OP, typename P1>
struct BiOperator<OP, P1, empty, empty> :
  public ParamVisitor<P1, SiconosVectorStorage&>
{
  typedef ParamVisitor<P1, SiconosVectorStorage&> base_type;

  BiOperator(const typename base_type::arguments_type& st) : base_type(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    apply_visitor<InternBiOp<OP, P1, const typename T::storage_type&> >
      (this->param2(), this->param1(), storage.internal_data);
  };
};



template<typename OP>
struct BiOperator<OP, empty, empty, empty> : public ParamVisitor<SiconosVectorStorage&>
{
  typedef ParamVisitor<SiconosVectorStorage&> base_type;

  /* the storage that is going to be modified is first kept in an argument */
  BiOperator(SiconosVectorStorage& st) : base_type(st) {};

  template<typename T>
  void operator() (const T& storage)
  {
    /* visit and modify param */
    apply_visitor<InternBiOp<OP, const typename T::storage_type&> >(this->param(), storage.internal_data);
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
  void operator()(const T2& x, const double& param, T1& y)
  {
    T2 tmp = param*x;
    y += tmp;
  }
};

struct OpAxpby
{
  template<typename T1, typename T2>
  void operator()(const T2& x, const double& a, const double& b, T1& y)
  {
    y *= b;
    y += a * x;
  }
};

template<>
void OpAxpby::operator()<DenseVect, DenseVect>(const DenseVect& x,
                                               const double& a,
                                               const double& b,
                                               DenseVect& y)
{
  bindings::scal(b, y);
  bindings::axpy(a, x, y);
}

template<>
void OpAxpby::operator()<SparseVect, SparseVect>(const SparseVect& x,
                                                 const double& a,
                                                 const double& b,
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
  void operator() (const T1& storage1, unsigned int param1, T2& storage2)
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
  void operator() (const T1& storage1, unsigned int param1, T2& storage2)
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
  void operator() (const T1& storage1, unsigned int param1, unsigned int param2, unsigned int param3, T2& storage2)
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
  void operator() (const T1& storage1, const double& param1, const Index& param2, bool param3, T2& storage2)
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
  Plus(SiconosVectorStorage& args) : base_type(args) {};
};
struct Minus : public BiOperator<OpMinus>
{
  typedef BiOperator<OpMinus> base_type;
  Minus(SiconosVectorStorage& args) : base_type(args) {};
};
struct Copy : public BiOperator<OpCopy>
{
  typedef BiOperator<OpCopy> base_type;
  Copy(SiconosVectorStorage& args) : base_type(args) {};
};
struct AddBlock : public BiOperator<OpAddBlock, unsigned int>
{
  typedef BiOperator<OpAddBlock, unsigned int> base_type;
  AddBlock(const base_type::arguments_type& args) : base_type(args) {};
};
struct SubBlock : public BiOperator<OpSubBlock, unsigned int>
{
  typedef BiOperator<OpSubBlock, unsigned int> base_type;
  SubBlock(const base_type::arguments_type& args) : base_type(args) {};
};
struct Axpy : public BiOperator<OpAxpy, double>
{
  typedef BiOperator<OpAxpy, double> base_type;
  Axpy(const base_type::arguments_type& args) : base_type(args) {};
};

struct Axpby : public BiOperator<OpAxpby, double, double>
{
  typedef BiOperator<OpAxpby, double, double> base_type;
  Axpby(const base_type::arguments_type& args) : base_type(args) {};
};

struct SetBlock : public BiOperator<OpSetBlock, unsigned int>
{
  typedef BiOperator<OpSetBlock, unsigned int> base_type;
  SetBlock(const base_type::arguments_type& args) : base_type(args) {};
};

struct Subscal : public BiOperator<OpSubscal, const double&, const Index&, const bool&>
{
  typedef BiOperator<OpSubscal, const double&, const Index&, const bool&> base_type;
  Subscal(const base_type::arguments_type& args) : base_type(args) {};
};
struct ToBlock : public BiOperator<OpToBlock, const unsigned int&, const unsigned int&, const unsigned int&>
{
  typedef BiOperator<OpToBlock, const unsigned int&, const unsigned int&, const unsigned int&> base_type;
  ToBlock(const base_type::arguments_type& args) : base_type(args) {};
};

struct StorageAllocator : public ParamVisitor<unsigned int>
{
  SiconosVectorStorage* answer;

  StorageAllocator(const unsigned int size) : ParamVisitor(size) {};

  template<typename T>
  void operator ()(const T&)
  {
    answer = new T(this->param());
  }
};


#endif
