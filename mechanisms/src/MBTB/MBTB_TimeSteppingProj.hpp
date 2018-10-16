#ifndef MBTB_TSPROJ_H
#define MBTB_TSPROJ_H
#include "SiconosKernel.hpp"
/**
 * \brief This class implements the time stepping with projection of a multi-bodies system.
 * It inherits from Siconos::TimeSteppingDirectProjection.
 * It consists in update the CAD word during the simulation.
 */
class MBTB_TimeSteppingProj : public TimeSteppingDirectProjection
{
public:
  /** Constructor with the time-discretisation.
   *  \param td pointer to a timeDiscretisation used in the integration
   *          (linked to the model that owns this simulation)
   *  \param osi one step integrator (default none)
   *  \param osnspb_velo one step non smooth problem (default none)
   *  \param osnspb_pos one step non smooth problem (default none)
   *  \param level
   */
  MBTB_TimeSteppingProj(
    SP::NonSmoothDynamicalSystem nsds,
    SP::TimeDiscretisation td,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb_velo,
    SP::OneStepNSProblem osnspb_pos,
    unsigned int level)
    :TimeSteppingDirectProjection(nsds, td,osi,osnspb_velo,osnspb_pos,level) {};
  //! Overloading of updateWorldFromDS.
  /*!
    It consists in updating the cad model from siconos.
   */
  virtual void updateWorldFromDS();  

};
TYPEDEF_SPTR(MBTB_TimeSteppingProj);
#endif
