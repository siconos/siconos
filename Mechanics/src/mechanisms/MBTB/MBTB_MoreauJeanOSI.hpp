#ifndef MBTB_MOREAU_H
#define MBTB_MOREAU_H
#include "SiconosKernel.hpp"
/**
 * \brief This class implements a variant of the std MoreauJeanOSI TS
 * It inherits from Siconos::MoreauJeanOSI
 * the main variants lies in the activation and desactivation of constraints
 */
class MBTB_MoreauJeanOSI : public MoreauJeanOSI
{
public:
  double _deactivateYPosThreshold;
  double _deactivateYVelThreshold;
  double _activateYPosThreshold;
  double _activateYVelThreshold;
public:
  /** constructor from a minimum set of data: one DS and its theta
   *  \param ds SP::DynamicalSystem : the DynamicalSystem linked to the OneStepIntegrator
   *  \param theta value
   */
  MBTB_MoreauJeanOSI(SP::DynamicalSystem ds , double theta);

  /** Apply the rule to one Interaction to known if is it should be included
   * in the IndexSet of level i
   * \param inter interaction
   * \param i level
   * \return a Boolean
   */
  bool addInteractionInIndexSet(SP::Interaction inter, unsigned int i);

  /** Apply the rule to one Interaction to known if is it should be removed
   * in the IndexSet of level i
   * \param inter interaction
   * \param i level
   * \return a Boolean
   */
  bool removeInteractionInIndexSet(SP::Interaction inter, unsigned int i);

};
TYPEDEF_SPTR(MBTB_MoreauJeanOSI);
#endif
