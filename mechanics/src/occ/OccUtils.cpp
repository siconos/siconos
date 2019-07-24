#include "OccUtils.hpp"
#include "OccContactFace.hpp"
#include "OccContactEdge.hpp"
#include "SiconosVector.hpp"
#include <TopoDS.hxx>
#include <gp_Dir.hxx>
#include <gp_Quaternion.hxx>
#include <BRepExtrema_DistShapeShape.hxx>

#include <cadmbtb.hpp>


void occ_move(TopoDS_Shape& shape, const SiconosVector& q)
{
  const gp_Vec translat = gp_Vec(q(0), q(1), q(2));
  const gp_Quaternion rota = gp_Quaternion(q(4), q(5), q(6), q(3));

  gp_Trsf transfo;
  transfo.SetRotation(rota);
  transfo.SetTranslationPart(translat);

  shape.Move(transfo);
  shape.Location(TopLoc_Location(transfo));
}

void occ_distanceFaceFace(const OccContactFace& csh1,
                          const OccContactFace& csh2,
                          Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                          Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                          Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                          Standard_Real& MinDist)
{
  // need the 2 sp pointers to keep memory
  SPC::TopoDS_Face pface1 = csh1.contact();
  SPC::TopoDS_Face pface2 = csh2.contact();

  const TopoDS_Face& face1 = *pface1;
  const TopoDS_Face& face2 = *pface2;

  BRepExtrema_DistShapeShape measure;
  measure.LoadS1(face1);
  measure.LoadS2(face2);
  measure.Perform();

  if (measure.IsDone())
  {
    /* we look for the first solution on a face */
    int nb_solutions = measure.NbSolution();
    for (Standard_Integer i=1; i<= nb_solutions; ++i)
    {
      if (measure.SupportTypeShape2(i) == BRepExtrema_IsInFace)
      {
         const gp_Pnt& p1 = measure.PointOnShape1(i);
         const gp_Pnt& p2 = measure.PointOnShape2(i);

         Standard_Real u, v;

         measure.ParOnFaceS2(i, u, v);
         gp_Dir normal = cadmbtb_FaceNormal(face2, u, v);
         normal.Coord(nX,nY,nZ);
         X1 = p1.X();
         X2 = p2.X();
         Y1 = p1.Y();
         Y2 = p2.Y();
         Z1 = p1.Z();
         Z2 = p2.Z();
         if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)<0)
         {
           normal.Reverse();
         }
         normal.Coord(nX,nY,nZ);
         MinDist = measure.Value();
         break;
      }
    }
  }
  else
    RuntimeException::selfThrow("occ distance: BRepExtrema_DistShapeShape failed");
}
void occ_distanceFaceEdge(const OccContactFace& csh1,
                          const OccContactEdge& csh2,
                          Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                          Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                          Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                          Standard_Real& MinDist)
{
  // need the 2 sp pointers to keep memory
  SPC::TopoDS_Face pface1 = csh1.contact();
  SPC::TopoDS_Edge pedge2 = csh2.contact();

  const TopoDS_Face& face1 = *pface1;
  const TopoDS_Edge& edge2 = *pedge2;

  BRepExtrema_DistShapeShape measure;
  measure.LoadS1(face1);
  measure.LoadS2(edge2);
  measure.Perform();

  if (measure.IsDone())
  {
    int nb_solutions = measure.NbSolution();
    for (Standard_Integer i=1; i<= nb_solutions; ++i)
    {
      /* we look for the first solution on a face */
      if (measure.SupportTypeShape1(i) == BRepExtrema_IsInFace)
      {
        const gp_Pnt& p1 = measure.PointOnShape1(i);
        const gp_Pnt& p2 = measure.PointOnShape2(i);

        Standard_Real u, v;

        measure.ParOnFaceS1(1, u, v);
        gp_Dir normal = cadmbtb_FaceNormal(face1, u, v);
        normal.Coord(nX,nY,nZ);
        X1 = p1.X();
        X2 = p2.X();
        Y1 = p1.Y();
        Y2 = p2.Y();
        Z1 = p1.Z();
        Z2 = p2.Z();
        if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)>0)
          normal.Reverse();
        normal.Coord(nX,nY,nZ);
        MinDist = measure.Value();
        break;
      }
    }
    // what to do now if MinDist is not changed ?
  }
  else
    RuntimeException::selfThrow("occ distance: BRepExtrema_DistShapeShape failed");
}
