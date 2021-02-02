#include <stdlib.h>
#include "lbl.h"

#ifdef __cplusplus
extern "C" {   /* To prevent C++ compilers from mangling symbols */
#endif

    // Safe wrapper around calloc().

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Calloc"
    void *LBL_Calloc(int nmemb, size_t size) {
        void *ptr;
        ptr = calloc(nmemb, size);
        if( ptr == NULL ) {
            SETERRQi(ALLOCATION_ERROR,
                     "Cannot allocate block of size", nmemb*size);
        }
        return ptr;
    }

    // Safe wrapper around realloc().

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Realloc"
    void *LBL_Realloc(void *ptr, int nmemb, size_t size) {
        ptr = realloc(ptr, nmemb*size);
        if( ptr == NULL ) {
            SETERRQi(REALLOCATION_ERROR,
                     "Cannot reallocate block of size", nmemb*size);
        }
        return ptr;
    }

    // Safe wrapper around free().
    // See also the LBL_Free shortcut in lbl.h.

#ifdef  __FUNCT__
#undef  __FUNCT__
#endif
#define __FUNCT__ "LBL_Free_Object"
    void LBL_Free_Object(void **ptr) {
        if( *ptr ) {
            free(*ptr);
            *ptr = NULL;
        }
        return;
    }

#ifdef __cplusplus
}
#endif
