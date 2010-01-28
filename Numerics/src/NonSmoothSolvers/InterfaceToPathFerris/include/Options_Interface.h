/*****************************************************************************/
/* Options_Interface.h                                                       */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Definitions needed to define an Options_Set.  Each option contains a    */
/*   name, a print field, a type, and a value.  An Options_Set contains an   */
/*   array of Options and callback routines.                                 */
/*                                                                           */
/*   NOTE: if you want to add your own options, contact the one of the       */
/*   authors of this code, and we will provide further information.          */
/*****************************************************************************/

#ifndef OPTIONS_INTERFACE_H
#define OPTIONS_INTERFACE_H

#include "Types.h"

/*****************************************************************************/
/* Constant declaration for the maximum size of the name and for the basic   */
/* types.                                                                    */
/*****************************************************************************/

#define OPTION_NAME_SIZE 128

#define Option_Boolean 0
#define Option_Integer 1
#define Option_Double 2

/*****************************************************************************/
/* Structure for the Option_Value.  NOTE: new option types from enumerated   */
/* data sets take integer values                                             */
/*****************************************************************************/

typedef struct
{
  Double d;
  Int i;
  Boolean b;
} Option_Value;

/*****************************************************************************/
/* Structure for the Option.                                                 */
/*****************************************************************************/

typedef struct
{
  Char name[OPTION_NAME_SIZE]; /* Name of the options                       */
  Boolean print;  /* Flag set to true if we want this options  */
  /* printed during a Options_Display( ).      */
  Int type;   /* Type for the option                       */
  Option_Value value;  /* Value for the option                      */
} Option;

/*****************************************************************************/
/* Structure for the Option_Set.                                             */
/*****************************************************************************/

struct _Option_Set
{
  Int numOptions;  /* Number of algorithm defined options       */
  Int numTypes;   /* Total number of types defined (3 plus any */
  /* new enumerated types)                     */

  Option *options;  /* Actual table of options                   */

  Void(CB_FPTR defaults)(Void); /* Function to set the defaults              */
  Void(CB_FPTR destroy)(Void); /* Function to destroy the option set        */

  Void(CB_FPTR *output)(const int, const int);
  /* Table for outputing the options.  First   */
  /* three elements are NULL, the rest are     */
  /* for the enumerated types you defined.     */
  /* First argument is the option number in    */
  /* the table, the second is the current      */
  /* value.                                    */

  Boolean(CB_FPTR *get_value)(const char *, int *);
  /* Table for obtaining the user defined      */
  /* options.  First three elements are NULL,  */
  /* the rest are for the enumerated types you */
  /* defined.  First argument is the string    */
  /* form of the option value, second argument */
  /* is the integer for enumerated type.       */
  /* Function returns True if the option value */
  /* was found in the table and false          */
  /* otherwise.                                */

  Void(CB_FPTR callback)(const int);
  /* Call back made when an option in this     */
  /* option set was successfully modified.     */

  Void(CB_FPTR common_callback)(const int);
  /* Call back made when an option in the      */
  /* common set of options was successfully    */
  /* modified.                                 */
};

#endif
