#ifndef MBTB_TSCOMBINEDPROJ_H
#define MBTB_TSCOMBINEDPROJ_H
#include "SiconosKernel.hpp"
/**
 * \brief This class implements the time stepping with projection of a multi-bodies system.
 * It inherits from Siconos::TimeSteppingProjectOnConstraints.
 * It consists in update the CAD word during the simulation.
 */
class MBTB_TimeSteppingCombinedProj : public TimeSteppingCombinedProjection
{
public:
  //! builder.
  MBTB_TimeSteppingCombinedProj(SP::TimeDiscretisation td,
                                SP::OneStepIntegrator osi,
                                SP::OneStepNSProblem osnspb_velo,
                                SP::OneStepNSProblem osnspb_pos,
                                unsigned int level) 
    :TimeSteppingCombinedProjection(td,osi,osnspb_velo,osnspb_pos,level) {} ;

  //! Overloading of updateWorldFromDS.
  /*!
    It consists in updating the cad model from siconos.
  */
  virtual void updateWorldFromDS();  
  
};
TYPEDEF_SPTR(MBTB_TimeSteppingCombinedProj);
#endif
