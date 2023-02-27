#include "Point.h"

#include <iostream>

Point::Point()
{
  x = 0.;
  y = 0.;
  z = 0.;
}

Point::Point(double a_x, double a_y, double a_z)
{
  x = a_x;
  y = a_y;
  z = a_z;
}

Point::~Point() {}

void Point::from_json(const json &j)
{
  j.at("x").get_to(x);
  j.at("y").get_to(y);
  j.at("z").get_to(z);
}

bool operator>(const Point &p1, const Point &p2) { return (p1.x >= p2.x); }

bool operator<(const Point &p1, const Point &p2) { return (p1.x < p2.x); }

void to_json(ojson &j, const Point &p)
{
  j.push_back(p.x);
  j.push_back(p.y);
  j.push_back(p.z);
}
