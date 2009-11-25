class SimpleMatrix;

class spSimpleMatrix : public boost::shared_ptr<SimpleMatrix>
{
public:
  spSimpleMatrix();
  spSimpleMatrix(SimpleMatrix * px): boost::shared_ptr<SimpleMatrix>(px)
  {
    ;
  }
  spSimpleMatrix(SimpleMatrix &x): boost::shared_ptr<SimpleMatrix>(&x, nullDeleter())
  {
    ;
  }
  spSimpleMatrix(int i);

};
class spConstSimpleMatrix : public boost::shared_ptr<const SimpleMatrix>
{
public:
  spConstSimpleMatrix();
  spConstSimpleMatrix(const SimpleMatrix &x): boost::shared_ptr<const SimpleMatrix>(&x, nullDeleter())
  {
    ;
  }
  spConstSimpleMatrix(int i);

};

typedef spSimpleMatrix  SPtrSimpleMatrix;
typedef spConstSimpleMatrix SPtrConstSimpleMatrix;

inline SPtrSimpleMatrix createSPtrSimpleMatrix(SimpleMatrix &x)
{
  spSimpleMatrix px(x);
  return px;
};
inline SPtrConstSimpleMatrix createSPtrConstSimpleMatrix(const SimpleMatrix &x)
{
  spConstSimpleMatrix px(x);
  return px;
} ;

namespace SharedPointer
{
typedef SPtrSimpleMatrix SimpleMatrix;
};

namespace SharedPointerConst
{
typedef SPtrConstSimpleMatrix X;
}
typedef ublas::compressed_matrix<SP::SimpleMatrix> BlocksSimpleMat;
