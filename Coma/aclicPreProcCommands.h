

#ifndef ACLICPREPROCCOMMANDS_H
#define ACLICPREPROCCOMMANDS_H


// ***********************************
//   Preprocessor Commands for ACLiC
// ***********************************
#define _OAWG
#define __STDC_LIMIT_MACROS
#define __STDC_CONSTANT_MACROS

// Define limits of integral types if using ACLiC
#if !defined(__CINT__) || defined(__MAKECINT__) 

#ifndef __APPLE__

/* Signed.  */
# define INT8_C(c)	c
# define INT16_C(c)	c
# define INT32_C(c)	c
# if __WORDSIZE == 64
#  define INT64_C(c)	c ## L
# else
#  define INT64_C(c)	c ## LL
# endif

/* Unsigned.  */
# define UINT8_C(c)	c
# define UINT16_C(c)	c
# define UINT32_C(c)	c ## U
# if __WORDSIZE == 64
#  define UINT64_C(c)	c ## UL
# else
#  define UINT64_C(c)	c ## ULL
# endif

/* Maximal type.  */
# if __WORDSIZE == 64
#  define INTMAX_C(c)	c ## L
#  define UINTMAX_C(c)	c ## UL
# else
#  define INTMAX_C(c)	c ## LL
#  define UINTMAX_C(c)	c ## ULL
# endif

#if __WORDSIZE == 64
#define __INT64_C(c)  c ## L
#define __UINT64_C(c) c ## UL
#else
#define __INT64_C(c)  c ## LL
#define __UINT64_C(c) c ## ULL
#endif                          //__WORDSIZE

/* Minimum of signed integral types.  */
#define INT8_MIN               (-128)
#define INT16_MIN              (-32767-1)
#define INT32_MIN              (-2147483647-1)
#define INT64_MIN              (-__INT64_C(9223372036854775807)-1)
/* Maximum of signed integral types.  */
#define INT8_MAX               (127)
#define INT16_MAX              (32767)
#define INT32_MAX              (2147483647)
#define INT64_MAX              (__INT64_C(9223372036854775807))
/* Maximum of unsigned integral types.  */
#define UINT8_MAX              (255)
#define UINT16_MAX             (65535)
#define UINT32_MAX             (4294967295U)
#define UINT64_MAX             (__UINT64_C(18446744073709551615))

#endif

// Have to hide VEGAS headers from CINT in interpreted mode
#include <VARootIO.h>
#include <VAShowerData.h>
#include <VAParameterData.h>
#endif      // !defined(__CINT__) || defined(_MAKECINT__)

#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <TH2F.h>
#include <TLine.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TPaveLabel.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TTimer.h>
#include <TPaveStats.h>
#include <TColor.h>
#endif
