#ifndef GEOMETER_HPP
#define GEOMETER_HPP

#include "MechanicsFwd.hpp"
#include "Question.hpp"
#include "ContactShapeDistance.hpp"
#include "OccUtils.hpp"
#include "cadmbtb.hpp"
#include <Standard_TypeDef.hxx>
#include <limits>
#include <iostream>

struct DistanceCalculatorType { VIRTUAL_ACCEPT_VISITORS(DistanceCalculatorType); };
struct OccDistanceType : DistanceCalculatorType { ACCEPT_STD_VISITORS(); };
struct CadmbtbDistanceType : DistanceCalculatorType { ACCEPT_STD_VISITORS(); };

struct Geometer : public Question<ContactShapeDistance>
{
  Geometer() {};
};

template<typename DistType>
void distanceFaceFace(const OccContactFace& csh1,
                      const OccContactFace& csh2,
                      Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                      Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                      Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                      Standard_Real& MinDist)
{}

template<typename DistType>
void distanceFaceEdge(const OccContactFace& csh1,
                      const OccContactEdge& csh2,
                      Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                      Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                      Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                      Standard_Real& MinDist)
{}

template<typename DistType>
void distanceEdgeEdge(const OccContactEdge& csh1,
                      const OccContactEdge& csh2,
                      Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                      Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                      Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                      Standard_Real& MinDist)
{
  throw "Geometer: Edge-Edge distance unimplemented";
}


template<>
void distanceFaceFace<CadmbtbDistanceType>(const OccContactFace& csh1,
                                   const OccContactFace& csh2,
                                   Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                   Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                   Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                   Standard_Real& MinDist)
{
  cadmbtb_distanceFaceFace(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ,
                           MinDist);
}

template<>
void distanceFaceEdge<CadmbtbDistanceType>(const OccContactFace& csh1,
                                   const OccContactEdge& csh2,
                                   Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                   Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                   Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                   Standard_Real& MinDist)
{
  cadmbtb_distanceFaceEdge(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ,
                           MinDist);
}


template<>
void distanceFaceFace<OccDistanceType>(const OccContactFace& csh1,
                                       const OccContactFace& csh2,
                                       Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                       Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                       Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                                       Standard_Real& MinDist)
{
  occ_distanceFaceFace(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ,
                       MinDist);
}

template<>
void distanceFaceEdge<OccDistanceType>(const OccContactFace& csh1,
                               const OccContactEdge& csh2,
                               Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                               Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                               Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                               Standard_Real& MinDist)
{
  occ_distanceFaceEdge(csh1, csh2, X1, Y1, Z1, X2, Y2, Z2, nX, nY, nZ, MinDist);
}

template <typename DistType>
struct FaceGeometer : public Geometer
{

  const OccContactFace& face1;

  FaceGeometer(const OccContactFace& face) : face1(face) {};
  using SiconosVisitor::visit;

  void visit(const OccContactFace& face2)
  {
    ContactShapeDistance& dist = this->answer;
    dist.value = std::numeric_limits<double>::infinity();
    distanceFaceFace<DistType>(this->face1, face2,
                               dist.x1, dist.y1, dist.z1,
                               dist.x2, dist.y2, dist.z2,
                               dist.nx, dist.ny, dist.nz,
                               dist.value);
  }
  void visit(const OccContactEdge& edge2)
  {
    ContactShapeDistance& dist = this->answer;
    dist.value = std::numeric_limits<double>::infinity();
    distanceFaceEdge<DistType>(this->face1, edge2,
                               dist.x1, dist.y1, dist.z1,
                               dist.x2, dist.y2, dist.z2,
                               dist.nx, dist.ny, dist.nz,
                               dist.value);
    dist.nx = -dist.nx;
    dist.ny = -dist.ny;
    dist.nz = -dist.nz;
  }

};

template<typename DistType>
struct EdgeGeometer : public Geometer
{

  const OccContactEdge& edge1;

  EdgeGeometer(const OccContactEdge& edge) : edge1(edge) {};
  using SiconosVisitor::visit;

  void visit(const OccContactFace& face2)
  {
    ContactShapeDistance& dist = this->answer;
    dist.value = std::numeric_limits<double>::infinity();
    distanceFaceEdge<DistType>(face2, this->edge1,
                               dist.x1, dist.y1, dist.z1,
                               dist.x2, dist.y2, dist.z2,
                               dist.nx, dist.ny, dist.nz,
                               dist.value);
  }
  void visit(const OccContactEdge& edge2)
  {
    ContactShapeDistance& dist = this->answer;
    dist.value = std::numeric_limits<double>::infinity();
    distanceEdgeEdge<DistType>(this->edge1, edge2,
                               dist.x1, dist.y1, dist.z1,
                               dist.x2, dist.y2, dist.z2,
                               dist.nx, dist.ny, dist.nz,
                               dist.value);
  }

};

#endif
