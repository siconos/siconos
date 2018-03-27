#ifndef OCC_UTILS
#define OCC_UTILS

#include "MechanicsFwd.hpp"
#include "SiconosFwd.hpp"
#include <Standard_TypeDef.hxx>

class OccContactFace;
class OccContactEdge;

void occ_move(TopoDS_Shape& shape, const SiconosVector& pos);

void occ_distanceFaceFace(const OccContactFace& csh1,
                          const OccContactFace& csh2,
                          Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                          Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                          Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                          Standard_Real& MinDist);

void occ_distanceFaceEdge(const OccContactFace& csh1,
                          const OccContactEdge& csh2,
                          Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                          Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                          Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                          Standard_Real& MinDist);

#endif
