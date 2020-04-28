
#ifndef MLCP_ENUM_TOOL_H
#define MLCP_ENUM_TOOL_H

/** Compute the total number of cases that should be enumerated
 * \param M the size of the MCLP problem.
 */
unsigned long long int computeNbCase(int M);


/** Initialize the enumeration process.
 * \param M the size of the MCLP problem.
 */
void initEnum(int M);

/** Iterate in the enumeration
 * \param[in,out] the next iterate
 */
int nextEnum(int * W2V);

#endif //MLCP_ENUM_H
