#ifndef SUPPORT_TYPEDEFINITIONS_H_
#define SUPPORT_TYPEDEFINITIONS_H_

/* *************************************
*  Compiler Specific Options
***************************************/
#ifdef _MSC_VER    /* Visual Studio */
#  pragma warning(disable : 4127)      /* disable: C4127: conditional expression is constant */
#  define FORCE_INLINE __forceinline
#else
#  if defined (__cplusplus) || defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L   /* C99 */
#    ifdef __GNUC__
#      define FORCE_INLINE inline __attribute__((always_inline))
#    else
#      define FORCE_INLINE inline
#    endif
#  else
#    define FORCE_INLINE inline
#  endif /* __STDC_VERSION__ */
#endif

//**************************************
// Basic Types
//**************************************
#if defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L   // C99
#include <stdint.h>
typedef int8_t		SBYTE;
typedef uint8_t		BYTE;
typedef int16_t		S16;
typedef uint16_t	U16;
typedef uint32_t	U32;
typedef int32_t		S32;
typedef uint64_t	U64;
typedef uint64_t	ULL;
typedef int64_t		S64;
#else
typedef char				SBYTE;
typedef unsigned char		BYTE;
typedef unsigned short		U16;
typedef short				S16;
typedef unsigned int		U32;
typedef signed int			S32;
typedef unsigned long long	U64;
typedef unsigned long long	ULL;
typedef signed long long	S64;
#endif

#endif /* SUPPORT_TYPEDEFINITIONS_H_ */
