/*! \addtogroup CADMBTB_PYTHON_API
  \brief This module contains the python API of the module CADMBTB.

  It provides the graphic functions allowing to plug an 3D view, from OCC. The 3D view is optional.<br>
  If there is no 3D view, the internal graphical management will be ignored.<br>
  It also provides a functions to load an artefact CAD model and set internal parameters. <br>
  For more details see the documentation of the functions.
  *  @{
  */


#ifndef CADMBTBPYTHONAPI
#define CADMBTBPYTHONAPI
/**
 * Load a step CAD file.
 * \param fileName const char * , a CAD file.
 * \param trans transparency value
 */
void CADMBTB_loadArtefactCADFile(const char * fileName, double trans);


class AIS_InteractiveContext;
class V3d_View;
/** It consists to initialize the graphic context.
 * (example : CADMBTB_setGraphicContext(display.Context) from PYTHON)
 * \param aisContext
 */
void CADMBTB_setGraphicContext(AIS_InteractiveContext & aisContext );

/** To set the current view.
 *  \param aView
 */
void CADMBTB_setGraphicView(V3d_View & aView);

/** To disable graphic 
 */
void CADMBTB_disableGraphic();

/** CADMBTB_setShapeDParam To set a double parameter.(extendable, without modifie the API)
 * This type of function has been chosen to easely set any parameters without modify the module API.
 * \param IdParam : identifier of the param.<br>
 *                0 for transparency.<br>
 * \param idShape : identifier of the shape.
 * \param v : value.
 */
void CADMBTB_setShapeDParam(unsigned int IdParam,unsigned int idShape,double  v);

/**To set a double parameter.(extendable, without modifie the API)
 *This type of function has been chosen to easily set any parameters without modifying 
 * the module API.
 * \param IdParam : identifier of the param.<br>
 *                0 for transparency.<br>
 * \param idContact : identifier of the contact.
 * \param idShape : identifier of the shape of the contact (0 or 1).
 * \param v : value.
 */
void CADMBTB_setContactDParam(unsigned int IdParam,unsigned int idContact,unsigned int idShape,double  v);

/** To set a int parameter. (extendable, without modifie the API)
 *  This type of function has been chosen to easely set any parameters without modifying
 * the module API.
 * \param IdParam : identifier of the param, used only for enable/disable the dump of the 3D view.
 * \param v : value.
 */
void CADMBTB_setIParam(unsigned int IdParam, int v);

/** Manual dump of the 3D view. */
void CADMBTB_DumpGraphic();


/** Print of the distance
 * \param v
 */

void CADMBTB_print_dist(unsigned int v);

/*! @} */
#endif
