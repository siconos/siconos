/* Siconos version 1.0, Copyright INRIA 2005.
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
#ifndef LAGRANGIANRELATION_H
#define LAGRANGIANRELATION_H

#include "Relation.h"
#include "LagrangianRXML.h"
#include "LagrangianDS.h"

/** \class LagrangianR
 *  \brief Lagrangian (Non Linear) Relation
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
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
 *   - holonom              -> h0(q), G0(q) (usual mechanical case)
 *   - non-holonom          -> h1(q,t), G10(q,t), G11(q,t)
 *   - holonom+control      -> h2(q,u), G20(q,u), G21(q,u)
 *   - non-holonom+control  -> h3(q,t,u), G30(q,t,u), G31(q,t,u), G32(q,t,u)
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
 * \warning At the time, only holonom and non-holonom constraints are modelled.
 *
 * Some rules:
 *   - h must be a plug-in
 *   - G can be either a fixed matrix or a plug-in, to fill this matrix in.
 *   - h and all the G functions must have the same signature.
 *   - G is defined as a stl-vector of SiconosMatrix* and the number of G functions (ie size of vector G)
 *     depends on the number of variables in h.
 *
 * Friction case: in Y, components are in the following order:
 *   first relation, normal part
 *   first relation, tangential part
 *   ...
 *   relation n, normal part
 *   relation n, tangential part
 * and so on ...
 * Note also that usually only normal part definition is required for Y[0].
 *
 */

class LagrangianRXML;

class LagrangianR : public Relation
{

protected:

  /* To define the type of constraints (holonom ...), ie the variables on which depend h and G*/
  std::string LagrangianRelationType;

  /* Boolean variables to check if h and G are plugged to external functions or not
   *  - always true for h
   *  - G can be either a plug-in or a given matrix
   *  - by default, h is not used in LagrangianLinear -> isHplugged = false
   *  - the same for G
   * That means that is..Plugged is false when connected to default plug-in.
   *  Note that these variables are mainly useful for derived classes.
   */
  bool isHPlugged;
  std::vector<bool> isGPlugged;

  /* Name of the plugin used to compute h */
  std::string  hFunctionName;

  /* G matrices that link Y[1] to dynamical system variables */
  std::vector<SiconosMatrix*> G;
  std::vector<bool> isGAllocatedIn;

  /* Name of the plug-in used to compute G */
  std::vector<std::string>  GFunctionName;

  // === plug-in, depending on problem type, ie LagrangianRelationType value ===

  // --- Plug-in for holonom case -> h0(q), G0(q) ---
  /** \fn void (*h0Ptr)(const unsigned int* sizeOfq, const double* q, const unsigned int* sizeOfy, const double* y);
   * \brief computes y = h(q)
   * \param unsigned int: sum of DS sizes, for DS involved in the interaction.
   * \param double*: vector q (composite of all DS q-vectors)
   * \param unsigned int: size of output vector y
   * \param double*: output vector y (in-out parameter)
   */
  void (*h0Ptr)(const unsigned int*, const double*, const unsigned int*, double*);
  /** \fn void (*G0Ptr)(const unsigned int* sizeOfq, const double* q, const unsigned int* sizeOfy, const double* G);
   * \brief computes G(q)
   * \param unsigned int: sum of DS sizes, for DS involved in the interaction.
   * \param double*: vector q (composite of all DS q-vectors)
   * \param unsigned int: size of output vector y
   * \param double*: output matrix G[0] (in-out parameter)
   */
  void (*G0Ptr)(const unsigned int*, const double*, const unsigned int*, double*);

  // --- plug-in for non-holonom case -> h1(q,t), G10(q,t), G11(q,t) ---
  /** \fn void (*h1Ptr)(const unsigned int* sizeOfq, const double* q,
   *                    const double* time, const unsigned int* sizeOfy, const double* y);
   * \brief computes y = h(q,t)
   * \param unsigned int: sum of DS sizes, for DS involved in the interaction.
   * \param double*: vector q (composite of all DS q-vectors)
   * \param double*: time
   * \param unsigned int: size of output vector y
   * \param double*: output vector y (in-out parameter)
   */
  void (*h1Ptr)(const unsigned int*, const double*, const double*, const unsigned int*, double*);
  // G0(q,t)
  /** \fn void (*G0Ptr)(const unsigned int* sizeOfq, const double* q,  const double* time,
   *                    const unsigned int* sizeOfy, const double* G);
   * \brief computes jacobian compare to q of function h
   * \param unsigned int: sum of DS sizes, for DS involved in the interaction.
   * \param double*: vector q (composite of all DS q-vectors)
   * \param double*: time
   * \param unsigned int: size of output vector y
   * \param double*: output matrix G[0] (in-out parameter)
   */
  void (*G10Ptr)(const unsigned int*, const double*, const double*, const unsigned int*, double*);
  // G1(q,t)
  /** \fn void (*G0Ptr)(const unsigned int* sizeOfq, const double* q, const double* time,
   *                    const unsigned int* sizeOfy, const double* G);
   * \brief computes jacobian compare to q of function h
   * \param unsigned int: sum of DS sizes, for DS involved in the interaction.
   * \param double*: vector q (composite of all DS q-vectors)
   * \param double*: time
   * \param unsigned int: size of output vector y
   * \param double*: output matrix G[1] (in-out parameter)
   */
  void (*G11Ptr)(const unsigned int*, const double*, const double*, const unsigned int*, double*);

  //\todo  ... other plug-in to be added if required ...

  /** class for plug-in management (open, close librairy...) */
  SiconosSharedLibrary cShared;

public:

  /** \fn LagrangianR(Interaction* =NULL);
   * \brief default constructor
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   */
  LagrangianR(Interaction* = NULL);

  /** \fn void LagrangianR(RelationXML*, Interaction* =NULL)
   *  \brief constructor from xml file
   *  \param relationXML
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  \exception RuntimeException
   */
  LagrangianR(RelationXML*, Interaction* = NULL);

  /** \fn void LagrangianR(const string& computeInput,const string& computeOutput, Interaction* =NULL)
   *  \brief constructor from a set of data
   *  \param string : the name of the plugin for computeInput
   *  \param string : the name of the plugin for computeOutput
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  \exception RuntimeException
   */
  LagrangianR(const std::string&, const std::string&, Interaction* = NULL);

  /** \fn LagrangianR(const Relation&)
   *  \brief copy constructor
   *  \param a relation to copy
   *  \param Interaction*: a pointer to the interaction that owns this relation (optional)
   *  warning: the interaction link is not copied, set a new one!
   */
  LagrangianR(const Relation &, Interaction* = NULL);

  ~LagrangianR();

  /** \fn  std::string getLagrangianRelationType() const
   *  \brief get the type of constraints of the relation (holonom ...)
   *  \return a string
   */
  inline const std::string getLagrangianRelationType() const
  {
    return LagrangianRelationType;
  }

  // -- G --

  /** \fn  vector<SiconosMatrix*> getGVector() const
   *  \brief get the vector of matrices G
   *  \return vector<SiconosMatrix*>
   */
  inline std::vector<SiconosMatrix*> getGVector() const
  {
    return G;
  }

  /** \fn void setGVector(const vector<SiconosMatrix*>& newVector)
   *  \brief set the value of G vector
   *  \param vector<SiconosMatrix*>
   */
  void setGVector(const std::vector<SiconosMatrix*> &);

  /** \fn const SiconosMatrix getG(unsigned int& index) const
   *  \brief get matrix G[index]
   *  \return a SiconosMatrix
   */
  inline const SiconosMatrix getG(const unsigned int & index = 0) const
  {
    return *(G[index]);
  }

  /** \fn SiconosMatrix* getGPtr() const
   *  \brief get a pointer on matrix G[index]
   *  \return a pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getGPtr(const unsigned int& index = 0) const
  {
    return G[index];
  }

  /** \fn void setG(const SiconosMatrix& newValue, unsigned int& index)
   *  \brief set the value of G[index] to newValue
   *  \param SiconosMatrix newValue
   *  \param unsigned int: index position in G vector
   */
  void setG(const SiconosMatrix&, const unsigned int& = 0);

  /** \fn void setGPtr(SiconosMatrix* newPtr,const unsigned int& index )
   *  \brief set G[index] to pointer newPtr
   *  \param SiconosMatrix * newPtr
   *  \param unsigned int: index position in G vector
   */
  void setGPtr(SiconosMatrix *newPtr, const unsigned int& = 0);


  /** \fn void setComputeHFunction(const string pluginPath, const string functionName&)
   *  \brief to set a specified function to compute function h(q,...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeHFunction(const std::string &, const std::string &);

  /** \fn void setComputeGFunction(const string pluginPath, const string functionName&, const unsigned int& index)
   *  \brief to set a specified function to compute G(q, ...)
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \param unsigned int: the index of G that must be computed (see introduction of this class for details on indexes)
   *  \exception SiconosSharedLibraryException
   */
  void setComputeGFunction(const std::string &, const std::string & , const unsigned int & = 0);

  /** \fn  std::string getHFunctionName() const
   *  \brief get name of function that computes h
   *  \return a string
   */
  inline const std::string getHFunctionName() const
  {
    return hFunctionName;
  }

  /** \fn  std::vector<string> getGFunctionNameVector() const
   *  \brief get list of names for functions that compute G
   *  \return a vector of strings
   */
  inline const std::vector<std::string> getGFunctionNameVector() const
  {
    return GFunctionName;
  }

  /** \fn  string getJacobianHFunctionName(const unsigned int& index) const
   *  \brief get names of functions that compute G[index]
   *  \return a string
   */
  inline std::string getGFunctionName(const unsigned int& index = 0) const
  {
    return GFunctionName[index];
  }

  /** \fn const bool isHPlugged() const
   *  \brief bool that shows if h is a plug-in or not
   *  \return a bool
   */
  //inline const bool isHPlugged() const {return isHPlugin;};

  /** \fn const vector<bool> isJacobianHPlugged() const
  *  \brief list of bool that shows if jacobianH is a plug-in or not
  *  \return a vector of bool
  */
  //inline const std::vector<bool> isJacobianHPlugged() const {return isJacobianHPlugin;};

  /** \fn void computeH(const double & time);
   * \brief to compute y = h(q,v,t) using plug-in mechanism
   * \param: double, current time
   */
  void computeH(const double &);

  /** \fn void computeG(const double & time, const unsigned int & index );
   * \brief to compute G using plug-in mechanism. Index shows which G is to be computed
   * \param: double, current time
   * \param: unsigned int
   */
  void computeG(const double &, const unsigned int& = 0);

  /** \fn void computeOutput(double time);
   *  \brief to compute output
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeOutput(const double&);

  /** \fn void computeFreeOutput(double time);
   *  \brief to compute y for the free state
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeFreeOutput(const double&);

  /** \fn void computeInput(double time);
   *  \brief to compute p
   *  \param double : current time
   *  \exception RuntimeException
   */
  virtual void computeInput(const double&);

  /** \fn void saveRelationToXML()
   *  \brief copy the data of the Relation to the XML tree
   */
  void saveRelationToXML() const;

  /** \fn LagrangianR* convert (Relation *r)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param Relation* : the relation which must be converted
   * \return a pointer on the relation if it is of the right type, NULL otherwise
   */
  static LagrangianR* convert(Relation *r);

  /** \fn  void display() const
   * \brief main relation members display
   */
  virtual void display() const;

};

#endif // LAGRANGIANRELATION_H
