//-----------------------------------------------------------------------------------------------------
// FILE:			CMassLynxRawScanReader.cpp
// AUTHOR:			Keith Richardson
// DATE:			Jan 2015
// COPYRIGHT(C):	Waters Corporation
//
// COMMENTS:        Pure "C" Wrapper around MassLynxRaw classes for use by MatLab / Mathematica 
//					
//-----------------------------------------------------------------------------------------------------
#include "MassLynxRawDefs.h"

//
// C Interface.
typedef void*   CMassLynxRawReader;
typedef void*   CMassLynxRawScanReader;

//
// Need an explicit constructor and destructor.
#ifdef __cplusplus
extern "C" {
#endif

// create a reader
MLYNX_RAW_API CMassLynxRawReader  newCMassLynxRawReader( 
    const char * cRawDataPath); // path to raw data

// delete a reader
MLYNX_RAW_API void   delCMassLynxRawReader(
    CMassLynxRawReader reader); // reader from newCMassLynxRawReader

// create a new scan reader
MLYNX_RAW_API CMassLynxRawScanReader  newCMassLynxRawScanReader(
    CMassLynxRawReader reader); // reader from newCMassLynxRawReader

// delete a scan reader
MLYNX_RAW_API void   delCMassLynxRawScanReader(
    CMassLynxRawScanReader scanreader); // scanreader from newCMassLynxRawScanReader

// get number of scans in a MassLynx function
MLYNX_RAW_API int getScansInFunction(
    CMassLynxRawScanReader scanreader, // scanreader from newCMassLynxRawScanReader
    int function );                    // MassLynx function number

// get number of points in a MassLynx scan
MLYNX_RAW_API int getScanSize(
    CMassLynxRawScanReader scanreader, // scanreader from newCMassLynxRawScanReader
    int function,                      // MassLynx function number
    int scan);                         // MassLynx scan number

// read a MassLynx spectrum
MLYNX_RAW_API int readSpectrum(
    CMassLynxRawScanReader scanreader, // scanreader from newCMassLynxRawScanReader 
    int function,                      // MassLynx function number
    int scan,                          // MassLynx scan number
    float * pMasses,                   // mass values (must be preallocated with correct length)
    float * pIntensities);             // intensities (must be preallocated with correct length)

// get X and y coordinates for specified function and scan number
MLYNX_RAW_API int getXYCoordinates( 
    CMassLynxRawReader reader, // reader from newCMassLynxRawReader
    int function,              // MassLynx function number
    int scan,                  // MassLynx scan number
    float * X,                 // pointer to X coordinate
    float * Y);                // pointer to Y coordinate

#ifdef __cplusplus
}
#endif
