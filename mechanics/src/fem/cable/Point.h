#pragma once
#include "BaseModel.h"
class Point : public BaseModel {
public:
  double x;
  double y;
  double z;

  Point();
  Point(double x, double y, double z);
  virtual ~Point();

  virtual void from_json(const json &j);

  Point &operator=(const Point &right)
  {
    x = right.x;
    y = right.y;
    z = right.z;
    return *this;
  }
  double norm() { return sqrt(x * x + y * y + z * z); }
  void diff(const Point &p1, const Point &p2)
  {
    x = p1.x - p2.x;
    y = p1.y - p2.y;
    z = p1.z - p2.z;
  }
  void add(const Point &p1, const Point &p2)
  {
    x = p1.x + p2.x;
    y = p1.y + p2.y;
    z = p1.z + p2.z;
  }
  void mult(double a_scalar)
  {
    x *= a_scalar;
    y *= a_scalar;
    z *= a_scalar;
  }
  void opposite()
  {
    x = -x;
    y = -y;
    z = -z;
  }
  void zero()
  {
    x = 0;
    y = 0;
    z = 0;
  }
  double dot(const Point &p) { return x * p.x + y * p.y + z * p.z; }

  friend bool operator>(const Point &p1, const Point &p2);
  friend bool operator<(const Point &p1, const Point &p2);
  friend void to_json(ojson &j, const Point &p);
};
