#ifndef _myExtrema_ExtSS_HeaderFile
#define _myExtrema_ExtSS_HeaderFile

#ifndef _Adaptor3d_SurfacePtr_HeaderFile
#include <Adaptor3d_SurfacePtr.hxx>
#endif
#ifndef _Standard_Boolean_HeaderFile
#include <Standard_Boolean.hxx>
#endif
#ifndef _Extrema_ExtElSS_HeaderFile
#include <Extrema_ExtElSS.hxx>
#endif
#ifndef _Extrema_SequenceOfPOnSurf_HeaderFile
#include <Extrema_SequenceOfPOnSurf.hxx>
#endif
#ifndef _Standard_Real_HeaderFile
#include <Standard_Real.hxx>
#endif
#ifndef _TColStd_SequenceOfReal_HeaderFile
#include <TColStd_SequenceOfReal.hxx>
#endif
#ifndef _GeomAbs_SurfaceType_HeaderFile
#include <GeomAbs_SurfaceType.hxx>
#endif
#ifndef _Standard_Integer_HeaderFile
#include <Standard_Integer.hxx>
#endif
class StdFail_NotDone;
class Standard_OutOfRange;
class Standard_TypeMismatch;
class Adaptor3d_Surface;
class Extrema_POnSurf;


#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif
/** 
 * \brief This class has been built from OCC in view of overloding the distance computation between CAD objects.
 */
//! It calculates all the extremum distances <br>
//!          between two surfaces. <br>
//!          These distances can be minimum or maximum. <br>
class myExtrema_ExtSS  {

public:

    void* operator new(size_t,void* anAddress) 
      {
        return anAddress;
      }
    void* operator new(size_t size) 
      { 
        return Standard::Allocate(size); 
      }
    void  operator delete(void *anAddress) 
      { 
        if (anAddress) Standard::Free((Standard_Address&)anAddress); 
      }
 // Methods PUBLIC
 // 


Standard_EXPORT myExtrema_ExtSS();

//! It calculates all the distances between S1 and S2. <br>
Standard_EXPORT myExtrema_ExtSS(const Adaptor3d_Surface& S1,const Adaptor3d_Surface& S2,const Standard_Real TolS1,const Standard_Real TolS2);

//! It calculates all the distances between S1 and S2. <br>
Standard_EXPORT myExtrema_ExtSS(const Adaptor3d_Surface& S1,const Adaptor3d_Surface& S2,const Standard_Real Uinf1,const Standard_Real Usup1,const Standard_Real Vinf1,const Standard_Real Vsup1,const Standard_Real Uinf2,const Standard_Real Usup2,const Standard_Real Vinf2,const Standard_Real Vsup2,const Standard_Real TolS1,const Standard_Real TolS2);

//! Initializes the fields of the algorithm. <br>
Standard_EXPORT   void Initialize(const Adaptor3d_Surface& S2,const Standard_Real Uinf2,const Standard_Real Usup2,const Standard_Real Vinf2,const Standard_Real Vsup2,const Standard_Real TolS1) ;

//! Computes the distances. <br>
//!          An exception is raised if the fieds have not been <br>
//!          initialized. <br>
Standard_EXPORT   void Perform(const Adaptor3d_Surface& S1,const Standard_Real Uinf1,const Standard_Real Usup1,const Standard_Real Vinf1,const Standard_Real Vsup1,const Standard_Real TolS1) ;

//! Returns True if the distances are found. <br>
Standard_EXPORT   Standard_Boolean IsDone() const;

//! Returns True if the curve is on a parallel surface. <br>
Standard_EXPORT   Standard_Boolean IsParallel() const;

//! Returns the number of extremum distances. <br>
Standard_EXPORT   Standard_Integer NbExt() const;

//! Returns the value of the Nth resulting distance. <br>
Standard_EXPORT   Standard_Real Value(const Standard_Integer N) const;

//! Returns the point of the Nth resulting distance. <br>
Standard_EXPORT   void Points(const Standard_Integer N,Extrema_POnSurf& P1,Extrema_POnSurf& P2) const;





protected:

 // Methods PROTECTED
 // 


 // Fields PROTECTED
 //


private: 

 // Methods PRIVATE
 // 


Standard_EXPORT   Adaptor3d_SurfacePtr Bidon() const;


 // Fields PRIVATE
 //
Adaptor3d_SurfacePtr myS2;
Standard_Boolean myDone;
Standard_Boolean myIsPar;
Extrema_ExtElSS myExtElSS;
Extrema_SequenceOfPOnSurf myPOnS1;
Extrema_SequenceOfPOnSurf myPOnS2;
Standard_Real myuinf1;
Standard_Real myusup1;
Standard_Real myvinf1;
Standard_Real myvsup1;
Standard_Real myuinf2;
Standard_Real myusup2;
Standard_Real myvinf2;
Standard_Real myvsup2;
Standard_Real mytolS1;
Standard_Real mytolS2;
TColStd_SequenceOfReal myvalue;
GeomAbs_SurfaceType myStype;


};





// other Inline functions and methods (like "C++: function call" methods)
//


#endif
