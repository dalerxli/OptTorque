/*
 * OptTorque.h -- a standalone code within the scuff-em suite 
 *                 -- for solving scattering problems
 *                 -- modeled after scuff-scatter.h 
 *
 * Homer Reid            -- 6/2011--2/2012
 * Yoonkyung Eunnie Lee  -- 2015.02.25
 */
#ifndef OPTTORQUE_H
#define OPTTORQUE_H 

#include "libhrutil.h"
#include "libhmat.h"
#include "libIncField.h"
#include "libscuff.h"

using namespace scuff;

/***************************************************************/
/* data structure containing everything needed to execute a    */
/* scattering calculation                                      */
/***************************************************************/
typedef struct SSData
 {
   RWGGeometry *G;
   HMatrix *M;
   HVector *RHS, *KN;
   cdouble Omega;
   double *kBloch;
   IncField *IF;
   double PowerRadius;
 } SSData;
 
/***************************************************************/
/* these are the 'output modules' that compute and process the */
/* scattered fields in various ways.                           */
/***************************************************************/
#if 0
void WriteOPFTFile(SSData *SSD, char *FileName, bool PlotFlux);
void WriteEPPFTFile(SSData *SSD, char *FileName, bool PlotFlux, int Order);
void WriteDSIPFTFile(SSData *SSD, char *FileName, char *DSIMesh,
                     double DSIRadius, int DSIPoints, bool DSICCQ,
                     bool DSIFarField, bool PlotFlux);
#endif
void WritePFTFile(SSData *SSD, PFTOptions *PFTOpts, int Method,
                  bool PlotFlux, char *FileName);
void WritePSDFile(SSData *SSD, char *PSDFile);
void GetMoments(SSData *SSD, char *MomentFile);
void ProcessEPFile(SSData *SSData, char *EPFileName);
void VisualizeFields(SSData *SSData, char *MeshFileName);
double GetIntegratedIntensity(RWGGeometry *G, int SurfaceIndex, HVector *RHSVector);

#endif
