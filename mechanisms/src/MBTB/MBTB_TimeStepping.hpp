#ifndef MBTB_TS_H
#define MBTB_TS_H
#include "SiconosKernel.hpp"
/**
 * \brief This class implements the time stepping of a multi-bodies system.
 * It inherits from Siconos::TimeStepping.
 * It consists in update the CAD word during the simulation.
 */
class MBTB_TimeStepping : public TimeStepping
{

public:
  /** Constructor with the time-discretisation.
   *  \param td pointer to a timeDiscretisation used in the integration
   *  (linked to the model that owns this simulation)
   *  \param osi one step integrator (default none)
   *  \param osnspb_velo one step non smooth problem (default none)
   */
  MBTB_TimeStepping(
    SP::NonSmoothDynamicalSystem nsds,
    SP::TimeDiscretisation td,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb_velo);
  //! Overloading of updateWorldFromDS.
  /*!
    It consists in updating the cad model from siconos.
   */
  virtual void updateWorldFromDS();

};
TYPEDEF_SPTR(MBTB_TimeStepping);
#endif
