/*****************************************************************************/
/* Import.h                                                                  */
/*                                                                           */
/* DESCRIPTION                                                               */
/*   Macro declarations for importing functions from a Windows DLL.  If you  */
/*   are using the windows DLL during your compilation, you MUST define      */
/*   USE_DLL.                                                                */
/*****************************************************************************/

#ifndef IMPORT_H
#define IMPORT_H

#if defined(_WIN32)
#  include <windows.h>
#  if defined(USE_EXPORT)
#    define FUN_DECL(type) __declspec(dllexport) type WINAPI
#  else
#    define FUN_DECL(type) type WINAPI
#  endif
#  define CB_FPTR WINAPI *
#  define CB_FUNC(type) type WINAPI
#  define SPD_SLEEP Sleep(1000)
#else
#  define FUN_DECL(type) type
#  define CB_FPTR *
#  define CB_FUNC(type) type
#  define SPD_SLEEP sleep(1)
#endif

#endif
