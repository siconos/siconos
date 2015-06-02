

#include "myExtrema_ExtSS.hxx"
#include "myExtrema_ExtSS.jxx"
#include "myExtrema_GenExtSS.hxx"
#include <StdFail_NotDone.hxx>
#include <Standard_NotImplemented.hxx>
#include <StdFail_InfiniteSolutions.hxx>
#include <Precision.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pln.hxx>

#include <Extrema_GenExtSS.hxx>
#include <ElCLib.hxx>

myExtrema_ExtSS::myExtrema_ExtSS() 
{
  myDone = Standard_False;
}

myExtrema_ExtSS::myExtrema_ExtSS(const Adaptor3d_Surface&     S1,
			     const Adaptor3d_Surface&     S2,
			     const Standard_Real    TolS1,
			     const Standard_Real    TolS2)

{
  Initialize(S2, S2.FirstUParameter(), 
	         S2.LastUParameter(), 
	         S2.FirstVParameter(), 
	         S2.LastVParameter(), TolS2);

  Perform(S1, S1.FirstUParameter(),
	  S1.LastUParameter(), 
	  S1.FirstVParameter(), 
	  S1.LastVParameter(), TolS1);
}

myExtrema_ExtSS::myExtrema_ExtSS(const Adaptor3d_Surface&   S1,
			     const Adaptor3d_Surface&   S2,
			     const Standard_Real  Uinf1,	
			     const Standard_Real  Usup1,
			     const Standard_Real  Vinf1,	
			     const Standard_Real  Vsup1,
			     const Standard_Real  Uinf2,	
			     const Standard_Real  Usup2,
			     const Standard_Real  Vinf2,	
			     const Standard_Real  Vsup2,
			     const Standard_Real  TolS1,
			     const Standard_Real  TolS2)

{
  Initialize(S2, Uinf2, Usup2, Vinf2, Vsup2, TolS2);
  Perform(S1, Uinf1, Usup1, Vinf1, Vsup1, TolS1);
}


void myExtrema_ExtSS::Initialize(const Adaptor3d_Surface&  S2,
			       const Standard_Real Uinf2,	
			       const Standard_Real Usup2,
			       const Standard_Real Vinf2,	
			       const Standard_Real Vsup2,
			       const Standard_Real TolS2)
{
  myS2 = (Adaptor3d_SurfacePtr)&S2;
  myIsPar = Standard_False;
  myuinf2  = Uinf2;
  myusup2  = Usup2;
  myvinf2  = Vinf2;
  myvsup2  = Vsup2;
  mytolS2  = TolS2;
  myStype  = S2.GetType();
}

				
void myExtrema_ExtSS::Perform(const Adaptor3d_Surface&   S1, 	
			    const Standard_Real  Uinf1,	
			    const Standard_Real  Usup1,
			    const Standard_Real  Vinf1,	
			    const Standard_Real  Vsup1,
			    const Standard_Real  TolS1)
{
  myuinf1  = Uinf1;
  myusup1  = Usup1;
  myvinf1  = Vinf1;
  myvsup1  = Vsup1;
  mytolS1 =  TolS1;
  myPOnS1.Clear();
  myPOnS2.Clear();
  myvalue.Clear();
  Standard_Integer i;
  //GeomAbs_SurfaceType myS1type  = S1.GetType();
  Standard_Integer NbU = 10, NbV = 10;
  
 
  myExtrema_GenExtSS Ext(S1, *myS2, NbU, NbV, mytolS1, mytolS2);
  myDone = Ext.IsDone();
  if (myDone) {
    Standard_Integer NbExt = Ext.NbExt();
    Standard_Real U1, V1,U2,V2;
    Extrema_POnSurf PS1;
    Extrema_POnSurf PS2;
    for (i = 1; i <= NbExt; i++) {
      PS1 = Ext.PointOnS1(i);
      PS2 = Ext.PointOnS2(i);
      PS1.Parameter(U1, V1);
      PS2.Parameter(U2, V2);
      if (S1.IsUPeriodic())
	U1 = ElCLib::InPeriod(U1, myuinf1, myuinf1+S1.UPeriod());
      if (S1.IsVPeriodic())
	V1 = ElCLib::InPeriod(V1, myvinf1, myvinf1+S1.VPeriod());
      if (myS2->IsUPeriodic())
	U2 = ElCLib::InPeriod(U2, myuinf2, myuinf2+myS2->UPeriod());
      if (myS2->IsVPeriodic())
	V2 = ElCLib::InPeriod(V2, myvinf2, myvinf2+myS2->VPeriod());
      
      if ((myuinf1-U1) <= mytolS1 && (U1-myusup1) <= mytolS1 &&
	  (myvinf1-V1) <= mytolS1 && (V1-myvsup1) <= mytolS1 &&
	  (myuinf2-U2) <= mytolS2 && (U2-myusup2) <= mytolS2 &&
	  (myvinf2-V2) <= mytolS2 && (V2-myvsup2) <= mytolS2) {
	myvalue.Append(Ext.Value(i));
	myPOnS1.Append(Extrema_POnSurf(U1, V1, PS1.Value()));
	myPOnS2.Append(Extrema_POnSurf(U2, V2, PS2.Value()));
      }
    }
  }
  
  
}


Standard_Boolean myExtrema_ExtSS::IsDone() const
{
  return myDone;
}

Standard_Boolean myExtrema_ExtSS::IsParallel() const
{
  return myIsPar;
}


Standard_Real myExtrema_ExtSS::Value(const Standard_Integer N) const
{
  if(!myDone) StdFail_NotDone::Raise();
  if (myIsPar && N != 1) StdFail_InfiniteSolutions::Raise();
  if ((N < 1) || (N > myvalue.Length())) Standard_OutOfRange::Raise();
  return myvalue.Value(N);
}


Standard_Integer myExtrema_ExtSS::NbExt() const
{
  if(!myDone) StdFail_NotDone::Raise();
  return myvalue.Length();
}



void myExtrema_ExtSS::Points(const Standard_Integer N,
			   Extrema_POnSurf&       P1,
			   Extrema_POnSurf&       P2) const
{
  if(!myDone) StdFail_NotDone::Raise();
  P1 = myPOnS1.Value(N);
  P2 = myPOnS2.Value(N);
}
