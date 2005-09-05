#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "NonSmoothDynamicalSystem.h"
#include "InteractionLink.h"
#include "Interaction.h"

// const
#include "SiconosConst.h"

#include <vector>
#include <string>
#include <map>

class NonSmoothDynamicalSystem;
class InteractionLink;
class Interaction;

/** \class Topology
 *  \brief this class provides maps to describe the topology of interactions of a NonSmoothDynamicalSystem
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) July 20, 2005
 *
 */

class Topology
{

private:

  // --- MEMBERS ---

  /** map that lists all the interactions and their linked interactions through common DS */
  std::map< Interaction*, std::vector<InteractionLink*>  > linkedInteractionMap;

  /** map that links interactions with their relative degrees */
  std::map< Interaction*, std::vector<unsigned int> > relativeDegreesMap;

  /** map that links interactions with the minimum index concerned by the nslaw in Y (the output derivatives vector of an interaction) */
  std::map< Interaction*, std::vector<unsigned int> > indexMinMap;

  /** map that links interactions with the maximum index concerned by the nslaw in Y (the output derivatives vector of an interaction)
   * indexMax <= relativeDegree - This map depends on the OneStepNSProblem*/
  std::map< Interaction*, std::vector<unsigned int> > indexMaxMap;

  /** global size of the effective output
   * effectiveSizeOutput = sum(indexMax - indexMin) over interactions */
  unsigned int effectiveSizeOutput;

  /** map that links each interaction with a list of indexes, giving the effective relations.
      An "effective" relation or output is one which is constrained.  */
  std::map< Interaction* , std::vector<unsigned int> > effectiveIndexesMap ;

  /** map that links interactions with their position in effective output vector */
  std::map< Interaction*, unsigned int> interactionEffectivePositionMap;

  /** check if topology has been updated since nsds modifications occur */
  bool isTopologyUpToDate;

  /** check if topology is static (all relative degrees = 0 or 1) or not */
  bool isTopologyTimeInvariant;

  /** the NonSmoothDynamicalSystem that owns this topology */
  NonSmoothDynamicalSystem * nsds;

  // === PRIVATE DEFAULT CONSTRUCTOR ===

  /** \fn Topology()
   *  \brief default constructor
   */
  Topology();

  // === OTHER PRIVATE FUNCTIONS ===

  /** \fn void computeLinkedInteractionMap()
   *   \brief compute the linkedInteractionMap
   * \param: a vector<Interaction*> (list of the interactions of the nsds)
   */
  void computeLinkedInteractionMap();

  /** \fn void computeRelativeDegreesMap()
   *   \brief compute the  RelativeDegreesMap
   */
  void computeRelativeDegreesMap();

  /** \fn vector<unsigned int> computeRelativeDegrees(Interaction *)
   *  \brief compute relative degrees vector of a specific interaction
   *  \param a pointer to Interaction
   */
  std::vector<unsigned int> computeRelativeDegrees(Interaction*);

  /** \fn void computeIndexMinMap()
   *   \brief compute the  IndexMinMap
   */
  void computeIndexMinMap();

  /** \fn vector<unsigned int> computeIndexMin(Interaction *)
   *  \brief compute relative degrees vector of a specific interaction
   *  \param a pointer to Interaction
   */
  std::vector<unsigned int> computeIndexMin(Interaction*);

public:

  // --- CONSTRUCTORS/DESTRUCTOR ---

  /** \fn Topology(NonSmoothDynamicalSystem*)
   *  \brief constructor from nsds that owns that topology
   * \param: a NonSmoothDynamicalSystem*
   */
  Topology(NonSmoothDynamicalSystem*) ;

  /** \fn Topology(const Topology&)
   *  \brief destructor */
  ~Topology();

  // === GETTERS/SETTERS ===

  // --- effectiveSizeOutput ---

  /** \fn const int getEffectiveSizeOutput() const
   *  \brief get the value of effectiveSizeOutput
   *  \return an unsigned int
   */
  inline const unsigned int getEffectiveSizeOutput() const
  {
    return effectiveSizeOutput;
  }

  /** \fn void setEffectiveSizeOutput(const int&)
   *  \brief set the value of effectiveSizeOutput
   *  \param an unsigned int
   */
  inline void setEffectiveSizeOutput(const unsigned int& newVal)
  {
    effectiveSizeOutput = newVal;
  }

  // --- linkedInteractionMap ---

  /** \fn  map< Interaction*, std::vector<InteractionLink*>> getLinkedInteractionMap(void)
   *  \brief get the linkedInteractionMap of this topology
   *  \return a map < Interaction*, std::vector<InteractionLink*>>
   */
  inline const std::map< Interaction*, std::vector<InteractionLink*> > getLinkedInteractionMap() const
  {
    return linkedInteractionMap;
  }

  // --- relativeDegreesMap ---

  /** \fn  map< Interaction*, std::vector<unsigned int> > getRelativeDegreesMap(void)
   *  \brief get the relativeDegreesMap of this topology
   *  \return a map < Interaction*, std::vector<unsigned int> >
   */
  inline const std::map< Interaction*, std::vector<unsigned int> > getRelativeDegreesMap() const
  {
    return relativeDegreesMap;
  }

  /** \fn  vector<unsigned int>  getRelativeDegrees(Interaction*)
   *  \brief get the relativeDegrees vector of a specific interaction
   *  \param a pointer on interaction
   *  \return a vector<unsigned int>
   */
  inline const std::vector<unsigned int> getRelativeDegrees(Interaction * Inter)
  {
    return relativeDegreesMap[Inter];
  }

  // --- indexMinMap ---

  /** \fn  map< Interaction*, std::vector<unsigned int> > getIndexMinMap(void)
   *  \brief get the indexMinMap of this topology
   *  \return a map < Interaction*, std::vector<unsigned int> >
   */
  inline const std::map< Interaction*, std::vector<unsigned int> > getIndexMinMap() const
  {
    return indexMinMap;
  }

  /** \fn  vector<unsigned int>  getIndexMin(Interaction*)
   *  \brief get the indexMin vector of a specific interaction
   *  \param a pointer on interaction
   *  \return a vector<unsigned int>
   */
  inline const std::vector<unsigned int> getIndexMin(Interaction * Inter)
  {
    return indexMinMap[Inter];
  }

  // --- indexMaxMap ---

  /** \fn  map< Interaction*, std::vector<unsigned int> > getIndexMaxMap(void)
   *  \brief get the indexMaxMap of this topology
   *  \return a map < Interaction*, std::vector<unsigned int> >
   */
  inline const std::map< Interaction*, std::vector<unsigned int> > getIndexMaxMap() const
  {
    return indexMaxMap;
  }

  /** \fn  vector<unsigned int>  getIndexMax(Interaction*)
   *  \brief get the indexMax vector of a specific interaction
   *  \param a pointer on interaction
   *  \return a vector<unsigned int>
   */
  inline const std::vector<unsigned int> getIndexMax(Interaction * Inter)
  {
    return indexMaxMap[Inter];
  }

  /** \fn  void  setIndexMax(Interaction*,vector<unsigned int> )
   *  \brief set the indexMax vector of a specific interaction
   *  \param a pointer on interaction
   *  \param a vector of int to set indexMax
   */
  inline void setIndexMax(Interaction * inter, const std::vector<unsigned int>&  index)
  {
    indexMaxMap[inter] = index;
  }

  // --- effectiveIndexesMap ---

  /** \fn  map< Interaction*, std::vector<unsigned int> > getEffectiveIndexesMap(void)
   *  \brief get the effectiveIndexesMap of this topology
   *  \return a map < Interaction*, std::vector<unsigned int> >
   */
  inline const std::map< Interaction*, std::vector<unsigned int> > getEffectiveIndexesMap() const
  {
    return effectiveIndexesMap;
  }

  /** \fn  vector<unsigned int>  getEffectiveIndexes(Interaction*)
   *  \brief get the effectiveIndexes vector of a specific interaction
   *  \param a pointer on interaction
   *  \return a vector<unsigned int>
   */
  inline const std::vector<unsigned int> getEffectiveIndexes(Interaction * Inter)
  {
    return effectiveIndexesMap[Inter];
  }

  /** \fn  void  setEffectiveIndexes(Interaction*,vector<unsigned int> )
   *  \brief set the effectiveIndexes vector of a specific interaction
   *  \param a pointer on interaction
   *  \param a vector of int to set effectiveIndexes
   */
  inline void setEffectiveIndexes(Interaction * inter, const std::vector<unsigned int>&  index)
  {
    effectiveIndexesMap[inter] = index;
  }

  // --- interactionEffectivePositionMap ---

  /** \fn  map<Interaction*, int> getOriginDSIndex(void)
   *  \brief get the interactionEffectivePositionMap of this topology
   *  \return a map <Interaction*, int>
   */
  inline const std::map< Interaction*, unsigned int> getInteractionEffectivePositionMap() const
  {
    return interactionEffectivePositionMap;
  }

  /** \fn  unsigned int getOriginDSIndex(void)
   *  \brief get the interactionEffectivePosition of a specific interaction
   *  \return an unsigned int
   */
  inline const unsigned int getInteractionEffectivePosition(Interaction * inter)
  {
    return interactionEffectivePositionMap[inter];
  }

  // --- isTopologyUpToDate ---

  /** \fn void  setUpToDate(const bool & val)
   *  \brief set isTopologyUpToDate to val
   *  \param a bool
   */
  inline void setUpToDate(const bool & val)
  {
    isTopologyUpToDate = val;
  }

  /** \fn bool isUpToDate()
   *  \brief check if topology has been updated since modifications occurs on nsds
   *  \return a bool
   */
  inline bool isUpToDate()
  {
    return isTopologyUpToDate;
  }

  // --- isTopologyTimeInvariant ---

  /** \fn void  setTimeInvariant(const bool & val)
   *  \brief set isTopologyTimeInvariant to val
   *  \param a bool
   */
  inline void setTimeInvariant(const bool & val)
  {
    isTopologyTimeInvariant = val;
  }

  /** \fn bool isTimeInvariant()
   *  \brief check if all relative degrees are equal to 0 or 1
   *  \return a bool
   */
  inline bool isTimeInvariant()
  {
    return isTopologyTimeInvariant;
  }

  /** \fn void updateTopology();
   *   \brief update topology: compute the linkedInteraction, position and relativeDegree maps and sizeOutput
   */
  void updateTopology();

  /** \fn void computeEffectiveSizeOutput()
   *   \brief compute effectiveSizeOutput, ie count the total number
   * of relations constrained
   */
  void computeEffectiveSizeOutput();

  /** \fn unsigned int computeEffectiveSizeOutput(Interaction *)
   *   \brief compute effectiveSizeOutput for a specific interaction
   * \return an unsigned int
   */
  unsigned int computeEffectiveSizeOutput(Interaction *);

  /** \fn void computeInteractionEffectivePositionMap()
   *   \brief compute the interactionEffectivePositionMap
   * \param: a vector<Interaction*> (list of the interactions of the nsds)
   */
  void computeInteractionEffectivePositionMap();

};

#endif // TOPOLOGY_H
