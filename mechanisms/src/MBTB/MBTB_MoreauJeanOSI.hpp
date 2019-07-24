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
   *  \param theta value for the theta parameter (default = 0.5)
   *  \param gamma value for the gamma parameter (default = NaN and gamma is not used)
   */
  MBTB_MoreauJeanOSI(double theta = 0.5 , double gamma = std::numeric_limits<double>::quiet_NaN());

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
  bool removeInteractionFromIndexSet(SP::Interaction inter, unsigned int i);

};
TYPEDEF_SPTR(MBTB_MoreauJeanOSI);
#endif
