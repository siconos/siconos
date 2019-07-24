/*! \addtogroup MBTB_INTERNAL_TOOL
   \brief This file contains the internal tools of the MBTB.

   * It consists in updating the CADMBTB from the simulation step. <br>
   * It also manages the output data.
   *  @{
   */
#ifndef INTERNALTOOLMBTB
#define INTERNALTOOLMBTB
#include <stdio.h>
#include <string>
#define PRINT_FORCE_CONTACTS
#define MBTB_PRINT_DIST
//! It updates the contacts CAD model from the body.
/*!

 */
void _MBTB_updateContactFromDS();

/** It updates the contacts CAD model from the body.
 * \param [in] numDS int,  update the cad model of contact related to the ds of id numDS.
 */
void _MBTB_updateContactFromDS(int numDS);

FILE* _MBTB_open(std::string filename, std::string args);

void _MBTB_close(FILE *);

/**!It prints the header of the output file.
 * \param fp output file
 */
void _MBTB_printHeader(FILE *fp);


/**It prints the current state in the output file.
 * \param fp output file
 */
void _MBTB_printStep(FILE *fp);

/** It displays the current state on std output.
 */
void _MBTB_displayStep();

/** It performs a step including the siconos call and the graphical update.
 */
void _MBTB_STEP();

#endif
/*! @} */
