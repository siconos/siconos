/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
/*! \file LagrangianR.h

*/
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.h"
#include "LagrangianRXML.h"
#include "LagrangianDS.h"


class LagrangianRXML;

//! Lagrangian (Non Linear) Relation
/** \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date Apr 27, 2004
 *  This class provides tools to describe non linear relation of the type:
 *
 * \f[
 * Y[0] = y = h(q,...)
 *
 * Y[1] = \dot y = G(q,...)\dot q + ...
 *
 * R = G^t(q,...)\lambda
 *
 * \f]
 *
 * using plug-in mechanism for functions h and G definition.
 *
 *
 * Recall:
 *   - output is vector Y, ie y and all its derivatives from 1 to r-1, r being the relative degree (see Interaction)
 *   - input is vector Lambda (see Interaction)
 *
 * member "LagrangianRelationType", set by user, defined h and G number of variables:
 *
 *   - scleronomic                     -> h0(q), G0(q) (usual mechanical case)
 *   - rhenomorous (ie time dependent) -> h1(q,t), G10(q,t), G11(q,t)
 *   - scleronomic+lambda              -> h2(q,lambda) , G20(q,lambda), G21(q,lambda)
 *   - scleronomic+control             -> h3(q,u), G30(q,u), G31(q,u)
 *   - rhenomourous+control            -> h4(q,t,u), G40(q,t,u), G41(q,t,u), G42(q,t,u)
 *
 * Usually, G functions corresponds to jacobians of h compare to its variables. But it is up to user
 * to define these plug-in.
 * Consider for example non-frictional case, with h(q,t).
 * Then:
 * \f[
 * Y[0] = y = h(q,t)
 *
 * Y[1] = \dot y = G[0](q,t)\dot q + G[1]
 *
 * with G[1] = dh/dt.
 *
 * \f]
 *
 *
 *
 * Depending on this type, computeH() function call the right plug-in. The same for G.
 * Mind that depending on the number of variables in h and G, corresponding jacobian operators
 * are to be defined if necessary.
 * \warning At the time, only scleronomic and non-scleronomic constraints are modelled.
 *
 * Some rules:
 *   - h must be a plug-in
 *   - G can be either a fixed matrix or a plug-in, to fill this matrix in.
 *   - h and all the G functions must have the same signature.
 *   - G is defined as a stl-vector of SiconosMatrix* and the number of G functions (ie size of vector G)
 *     depends on the number of variables in h.
 *
 * Friction case: in Y, components are in the following order:\n
 *   first relation, normal part\n
 *   first relation, tangential part\n
 *   ...\n
 *   relation n, normal part \n
 *   relation n, tangential part\n
 * and so on ...\n
 * Note also that usually only normal part definition is required for Y[0].
 *
 */
class LagrangianR : public Relation
{

protected:

  /** To define the type of constraints (scleronomic ...), ie the variables on which depend h and G*/
  std::string LagrangianRelationType;

  /** Boolean variables to check if h and G are plugged to external functions or not
   *  - always true for h
   *  - G can be either a plug-in or a given matrix
   *  - by default, h is not used in LagrangianLinear -> isHplugged = false
   *  - the same for G
   * That means that is..Plugged is false when connected to default plug-in.
   *  Note that these variables are mainly useful for derived classes.
   */
  bool isHPlugged;
  std::deque<bool> isGPlugged;

  /** Name of the plugin used to compute h */
  std::string  hFunctionName;

  /** G matrices that link Y[1] to dynamical system variables */
  std::vector<SiconosMatrix*> G;

  /** Name of the plug-in used to compute G */
  std::vector<std::string>  GFunctionName;

  // === plug-in, depending on problem type, ie LagrangianRelationType value ===

  // --- Plug-in for scleronomic case -> h0(q), G0(q) ---
  /** LagrangianR plug-in to compute h0(q) (scleronomic case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param sizeY : size of vector y (ie of the intercation)
   * @param[in,out] y : pointer to the first element of y
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*h0Ptr)(unsigned int, const double*, unsigned int, double*, double*);

  /** LagrangianR plug-in to compute G0(q), gradient of h0 according to q (scleronomic case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param sizeY : size of vector y (ie of the intercation)
   * @param[in,out] G0 : pointer to the first element of G0 (sizeY X sizeDS matrix)
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*G0Ptr)(unsigned int, const double*, unsigned int, double*, double*);

  // --- plug-in for non-scleronomic case -> h1(q,t), G10(q,t), G11(q,t) ---
  /** LagrangianR plug-in to compute h1(q,t) (non-scleronomic case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param time : current time
   * @param sizeY : size of vector y (ie of the intercation)
   * @param[in,out] y : pointer to the first element of y
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*h1Ptr)(unsigned int, const double*, double, unsigned int, double*, double*);

  /** LagrangianR plug-in to compute G10(q,t), gradient of h1 accoring to q (non-scleronomic case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param time : current time
   * @param sizeY : size of vector y (ie of the intercation)
   * @param[in,out] G10 : pointer to the first element of G10
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*G10Ptr)(unsigned int, const double*, double, unsigned int, double*, double*);

  /** LagrangianR plug-in to compute G11(q,t), gradient of h1 according to time (non-scleronomic case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param time : current time
   * @param sizeY : size of vector y (ie of the intercation)
   * @param[in,out] G11 : pointer to the first element of G11
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*G11Ptr)(unsigned int, const double*, double, unsigned int, double*, double*);

  // --- plug-in for scleronomic+lambda case -> h2(q,lambda), G20(q,lambda), G21(q,lambda) ---

  /** LagrangianR plug-in to compute h2(q,lambda) (scleronomic+lambda case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param lambda : current time
   * @param sizeY : size of vector y (ie of the interaction)
   * @param[in,out] y : pointer to the first element of y
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*h2Ptr)(unsigned int, const double*, const double*, unsigned int, double*, double*);

  /** LagrangianR plug-in to compute G20(q,lambda), gradient of h2 according to q (scleronomic+lambda case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param lambda : current time
   * @param sizeY : size of vector y (ie of the interaction)
   * @param[in,out] G20 : pointer to the first element of G20
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*G20Ptr)(unsigned int, const double*, const double*, unsigned int, double*, double*);

  /** LagrangianR plug-in to compute G21(q,lambda), gradient of h2 according to lambda (scleronomic+lambda case)
   * @param sizeDS : sum of the sizes of all the DynamicalSystems involved in the interaction
   * @param q : pointer to the first element of q
   * @param lambda : current time
   * @param sizeY : size of vector y (ie of the interaction)
   * @param[in,out] G21 : pointer to the first element of G21
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*G21Ptr)(unsigned int, const double*, const double*, unsigned int, double*, double*);

  //\todo  ... other plug-in to be added if required ...

  /** class for plug-in management (open, close librairy...) */
  SiconosSharedLibrary cShared;

public:

  /** default constructor
   */
  LagrangianR();

  /** constructor from xml file
   *  \param relationXML
   *  \exception RuntimeException
   */
  LagrangianR(RelationXML*);

  /** constructor from a set of data
   *  \param string : the type of relation (scleronomic ...)
   *  \param string : the name of the plugin for computeH
   *  \param vector<string> : a list of names for the plugin for computeG (depends on relation type)
   *  \exception RuntimeException
   */
  LagrangianR(const std::string, const std::string, const std::vector<std::string>&);

  /** copy constructor
   *  \param a relation to copy
   *  warning: the interaction link is not copied, set a new one!
   */
  LagrangianR(const Relation &);

  /** destructor
   */
  virtual ~LagrangianR();

  /** initialize the relation (check sizes, memory allocation ...)
   */
  virtual void initialize();

  /** get the type of constraints of the relation (scleronomic ...)
   *  \return a string
   */
  inline const std::string getLagrangianRelationType() const
  {
    return LagrangianRelationType;
  }

  /** check and/or allocate memory for G
   */
  void manageGMemory();

  /** set the type of constraints of the relation (scleronomic ...) and adapt corresponding variables
   * (resize G ...)
   *  \param a string
   */
  void setLagrangianRelationType(const std::string);

  // -- G --

  /** get the vector of matrices G
   *  \return vector<SiconosMatrix*>
   */
  inline std::vector<SiconosMatrix*> getGVector() const
  {
    return G;
  }

  /** set the value of G vector
   *  \param vector<SiconosMatrix*>
   */
  void setGVector(const std::vector<SiconosMatrix*> &);

  /** get matrix G[index]
   *  \return a SimpleMatrix
   */
  inline const SimpleMatrix getG(const unsigned int  index = 0) const
  {
    return *(G[index]);
  }

  /** get a pointer on matrix G[index]
   *  \return a pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getGPtr(const unsigned int index = 0) const
  {
    return G[index];
  }

  /** set the value of G[index] to newValue
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in G vector
   */
  void setG(const SiconosMatrix&, const unsigned int = 0);

  /** set G[index] to pointer newPtr
   *  \param SiconosMatrix * newPtr
   *  \param unsigned int: index position in G vector
   */
  void setGPtr(SiconosMatrix *newPtr, const unsigned int = 0);


  /** to set a specified function to compute function h(q,...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeHFunction(const std::string , const std::string);

  /** to set a specified function to compute G(q, ...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \param unsigned int: the index of G that must be computed (see introduction of this class for details on indexes)
   *  \exception SiconosSharedLibraryException
   */
  void setComputeGFunction(const std::string , const std::string  , const unsigned int  = 0);

  /** get name of function that computes h
   *  \return a string
   */
  inline const std::string getHFunctionName() const
  {
    return hFunctionName;
  }

  /** get list of names for functions that compute G
   *  \return a vector of strings
   */
  inline const std::vector<std::string> getGFunctionNameVector() const
  {
    return GFunctionName;
  }

  /** get names of functions that compute G[index]
   *  \return a string
   */
  inline std::string getGFunctionName(const unsigned int index = 0) const
  {
    return GFunctionName[index];
  }

  /** to compute y = h(q,v,t) using plug-in mechanism
   * \param: double, current time
   */
  void computeH(const double);

  /** to compute G using plug-in mechanism. Index shows which G is to be computed
   * \param: double, current time
   * \param: unsigned int
   */
  void computeG(const double , const unsigned int = 0);

  /** to compute output
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeOutput(const double, const unsigned int = 0);

  /** to compute output y[0]
   *  \param double : current time
   *  \param a pointer to SiconosVector: q for all the DS
   */
  void computeY0(const double, SiconosVector*);

  /** to compute output y[1]
   *  \param double : current time
   *  \param a pointer to SiconosVector: q for all the DS
   *  \param a pointer to SiconosVector: velocity for all the DS
   */
  void computeY1(const double, SiconosVector*, SiconosVector*);

  /** to compute output y[2]
   *  \param double : current time
   *  \param a pointer to SiconosVector: q for all the DS
   *  \param a pointer to SiconosVector: acceleration for all the DS
   */
  void computeY2(const double, SiconosVector*, SiconosVector*);

  /** to compute y for the free state
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute, optional, default = 0.
   */
  virtual void computeFreeOutput(const double, const unsigned int = 0);

  /** to compute p
   *  \param double : current time
   *  \param unsigned int: "derivative" order of lambda used to compute input
   */
  virtual void computeInput(const double, const unsigned int);

  /** get in Matrix G[index] the block corresponding to DS number int
   *  \param int, the ds number
   *  \param SiconosMatrix (in-out parameter): the resulting block matrix
   */
  void getGBlockDS(const int, SiconosMatrix&, const unsigned int = 0) const;


  /** get in Matrix G[index] the block corresponding to ds
   *  \param a pointer to a dynamical system
   *  \param SiconosMatrix (in-out parameter): the resulting block matrix
   */
  void getGBlockDS(DynamicalSystem *, SiconosMatrix&, const unsigned int = 0) const;

  /** copy the data of the Relation to the XML tree
   */
  void saveRelationToXML() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation* : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianR* convert(Relation *r);

  /** main relation members display
   */
  virtual void display() const;

};

#endif // LAGRANGIANRELATION_H
