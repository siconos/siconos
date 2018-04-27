#ifndef cadmbtb_hpp
#define cadmbtb_hpp

class gp_Pnt;
class gp_Dir;
class TopoDS_Face;
class TopoDS_Edge;

#include <MechanicsFwd.hpp>
#include <Standard_TypeDef.hxx>

gp_Pnt cadmbtb_FacePoint(const TopoDS_Face &face,Standard_Real u, Standard_Real v);
gp_Pnt cadmbtb_EdgePoint(const TopoDS_Edge &edge,Standard_Real u);
gp_Dir cadmbtb_FaceNormal(const TopoDS_Face &face,Standard_Real u, Standard_Real v);

void cadmbtb_distanceFaceFace(const OccContactFace& csh1,
                              const OccContactFace& csh2,
                              Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                              Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                              Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                              Standard_Real& MinDist);

void cadmbtb_distanceFaceEdge(
  const OccContactFace& sh1, const OccContactEdge& sh2,
  Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
  Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
  Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
  Standard_Real& MinDist);

#endif
