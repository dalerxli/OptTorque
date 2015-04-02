/***************************************************************
 * OptTorque.cc -- Compute Power,Force,Torque 
 * Last updated on 2015.03.15, v11
 *
 ***************************************************************
 * v5 & v55 & v6: Functional version of OptTorque.cc 
 *                but only GHBeam/GLBeam sources were allowed 
 * updates to v7: Implemented all of scuff-scatter functionalities
 * updates to v8: GHBeam Declaration with correct pointers
 * updates to v9: SSD->IF=IFDList; Multiple GBeam Sources allowed 
 * updates to v10: Gouy Phase corrected
 * updates to v11: GetIntegratedIntensity added 
 * updates to v12: Intensity multiplier is added 
 *
 ***************************************************************/
#include <stdio.h>
#include <math.h>
#include <complex>
#include <stdarg.h>
#include <fenv.h>

#include "GHBeam.h"
#include "OptTorque.h"

#define II cdouble(0.0,1.0)
/***************************************************************/
/***************************************************************/
#define MAXPW    10    // max number of plane waves
#define MAXGB    1     // max number of gaussian beams
#define MAXGLB   10    // max number of GLBeams 
#define MAXGHB   10    // max number of GHBeams
#define MAXPS    10    // max number of point sources
#define MAXFREQ  10    // max number of frequencies
#define MAXEPF   10    // max number of evaluation-point files
#define MAXFVM   10    // max number of field visualization meshes
#define MAXCACHE 10    // max number of cache files for preload

#define MAXSTR   1000

using namespace scuff;

/***************************************************************/
/***************************************************************/
/***************************************************************/
// These static declarations are redundant. 
static char *FieldFuncs=const_cast<char *>(
 "|Ex|,|Ey|,|Ez|,"
 "sqrt(|Ex|^2+|Ey|^2+|Ez|^2),"
 "|Hx|,|Hy|,|Hz|,"
 "sqrt(|Hx|^2+|Hy|^2+|Hz|^2)");

static const char *FieldTitles[]=
 {"|Ex|", "|Ey|", "|Ez|", "|E|",
  "|Hx|", "|Hy|", "|Hz|", "|H|",
 };

#define NUMFIELDFUNCS 8
/***************************************************************/
/***************************************************************/
/***************************************************************/
//void VisualizeFields(RWGGeometry *G, IncField *IF, HVector *KN,
//                     cdouble Omega, char *MeshFileName);
//double **AllocateByEdgeArray(RWGGeometry *G, int ns);
//void ProcessByEdgeArray(RWGGeometry *G, int ns, cdouble Omega,
//                        double **ByEdge);

/***************************************************************/
/***************************************************************/
/***************************************************************/
int main(int argc, char *argv[])
{  
  /***************************************************************/
  /* process options *********************************************/
  /***************************************************************/
  InstallHRSignalHandler();
//
  char *GeoFile=0;
//
  cdouble OmegaVals[MAXFREQ];        int nOmegaVals;
  char *OmegaFile=0;
//
  double pwDir[3*MAXPW];             int npwDir;
  cdouble pwPol[3*MAXPW];            int npwPol;
//
  double gbDir[3*MAXGB];             int ngbDir;
  cdouble gbPol[3*MAXGB];            int ngbPol;
  double gbCenter[3*MAXGB];          int ngbCenter;
  double gbWaist[MAXGB];             int ngbWaist;
//
  double ghbWaist[MAXGHB];            int nghbWaist; 
  int HM[MAXGHB];                     int ngHM;
  int HN[MAXGHB];                     int ngHN; 

  double glbWaist[MAXGLB];            int nglbWaist;
  char *glbpolname[MAXGLB];              int nglbpolname; 
  double glbCenter[3*MAXGLB];          int nglbCenter;
  double glbDir[3*MAXGLB];             int nglbDir;
  int P[MAXGLB];                      int nglP; 
  int L[MAXGLB];                      int nglL; 
  double glbI0[MAXGLB];               int nglbI0; 
//
  double psLoc[3*MAXPS];             int npsLoc;
  cdouble psStrength[3*MAXPS];       int npsStrength;
//
  char *FVMeshes[MAXFVM];            int nFVMeshes;
//
  char *EPFiles[MAXEPF];             int nEPFiles;
//
  char *PFTFile=0;
  char *OPFTFile=0;
  char *EPPFTFile=0;
  int  EPFTOrder=1;
//
  char *DSIPFTFile = 0;
  double DSIRadius = 10.0;
  int DSIPoints    = 302;
  char *DSIMesh    = 0;
  bool DSIFarField = false;
//
  bool PlotPFTFlux=false;
//
  char *MomentFile=0;
  char *PSDFile=0;
  bool PlotSurfaceCurrents=false;
//
  bool MatrixOut=false; ///
  char *HDF5File=0;
  char *Cache=0;
  char *ReadCache[MAXCACHE];         int nReadCache;
  char *WriteCache=0;
  char *LogLevel=0;
//
  char *FileBase=0;
//
  /* name        type    #args  max_instances  storage    count  description*/
  OptStruct OSArray[]=
   {
     {"geometry",  PA_STRING,  1, 1, (void *)&GeoFile,   0,  ".scuffgeo file"},
     {"Omega",     PA_CDOUBLE, 1, MAXFREQ, (void *)OmegaVals, &nOmegaVals,  "(angular) frequency"},
     {"OmegaFile", PA_STRING,  1, 1, (void *)&OmegaFile, 0,  "file listing angular frequencies"},
     {"FileBase",  PA_STRING,  1, 1, (void *)&FileBase,      0,  "base filename for output file"},
//
     {"gbDirection",PA_DOUBLE, 3, MAXGB, (void *)gbDir,     &ngbDir,   "gaussian beam direction"},
     {"gbPolarization", PA_CDOUBLE, 3, MAXGB, (void *)gbPol,&ngbPol,   "gaussian beam polarization"},
     {"gbCenter",  PA_DOUBLE,  3, MAXGB, (void *)gbCenter,  &ngbCenter,"gaussian beam center"},
     {"gbWaist",   PA_DOUBLE,  1, MAXGB, (void *)gbWaist,   &ngbWaist, "gaussian beam waist"},
//
     {"pwDirection",  PA_DOUBLE,  3, MAXPW, (void *)pwDir,   &npwDir,  "plane wave direction"},
     {"pwPolarization", PA_CDOUBLE, 3, MAXPW, (void *)pwPol, &npwPol,  "plane wave polarization"},
//
//     {"GBeamMode", PA_INT,     1, 1, (void *)&GBeamMode,     0,  "Beam Mode Index, 1:Gauss-Hermite, 2:Gauss-Legendre"},
//
     {"HM",        PA_INT,     1, MAXGHB, (void *)HM,       &ngHM,  "Hermite(x) M parameter"},
     {"HN",        PA_INT,     1, MAXGHB, (void *)HN,       &ngHN,  "Hermite(y) N parameter"},
     {"ghbWaist",  PA_DOUBLE,  1, MAXGHB, (void *)ghbWaist, &nghbWaist, "Hermite beam waist"},
//
     {"P",         PA_INT,     1, MAXGLB, (void *)P,        &nglP,  "Laguerre P(radial) parameter"},
     {"L",         PA_INT,     1, MAXGLB, (void *)L,        &nglL,  "Laguerre L(azimuthal) parameter"},
     {"glbWaist",  PA_DOUBLE,  1, MAXGLB, (void *)glbWaist, &nglbWaist, "Laguerre beam waist"},
     {"glbI0",     PA_DOUBLE,  1, MAXGLB, (void *)glbI0, &nglbI0, "Laguerre beam Intensity Factor"},
     {"glbpolname", PA_STRING, 1, MAXGLB, (void *)&glbpolname, &nglbpolname, "Laguerre beam polarization name, one of : xp, yp, rcp, lcp, radial, azimuthal "},
     {"glbCenter",  PA_DOUBLE,  3, MAXGLB, (void *)glbCenter,  &nglbCenter,"Laguerre beam center"},
     {"glbDirection",PA_DOUBLE, 3, MAXGLB, (void *)glbDir,&nglbDir,"Laguerre beam direction"},
//
     {"psLocation",PA_DOUBLE,  3, MAXPS,(void *)psLoc,     &npsLoc,       "point source location"},
     {"psStrength",PA_CDOUBLE, 3, MAXPS,(void *)psStrength,&npsStrength,  "point source strength"},
//
     {"MatrixOut", PA_BOOL,    0, 1, (void *)&MatrixOut,     0,  "output F, M, KN into .dat files"},
     {"PlotPFTFlux",PA_BOOL,    0, 1, (void *)&PlotPFTFlux,      0,  "generate plots of spatially-resolved PFT flux"},
     {"EPFile",    PA_STRING,  1, MAXEPF,(void *)EPFiles, &nEPFiles,  "list of evaluation points"},
//
     {"FVMesh",    PA_STRING,  1, MAXFVM,(void *)FVMeshes,&nFVMeshes,  "field visualization mesh"},

     {"PFTFile",   PA_STRING,  1, 1, (void *)&PFTFile,    0, "name of PFT output file"},
     {"OPFTFile",  PA_STRING,  1, 1, (void *)&OPFTFile,   0, "name of overlap PFT output file"},
     {"EPPFTFile", PA_STRING,  1, 1, (void *)&EPPFTFile,  0, "name of equivalence-principle PFT output file"},
     {"EPFTOrder", PA_INT,     1, 1, (void *)&EPFTOrder,  0, "cubature order for equivalence-principle force/torque (1,4,9,13,20)"},
     {"DSIPFTFile",PA_STRING,  1, 1, (void *)&DSIPFTFile, 0, "name of displaced surface-integral PFT output file"},
     {"DSIMesh",   PA_STRING,  1, 1, (void *)&DSIMesh,    0, "mesh file for surface-integral PFT"},
     {"DSIRadius", PA_DOUBLE,  1, 1, (void *)&DSIRadius,  0, "radius of bounding sphere for DSIPFT"},
     {"DSIPoints", PA_INT,     1, 1, (void *)&DSIPoints,  0, "number of quadrature points for surface-integral PFT (6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434, 590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810)"},
     {"DSIFarField", PA_BOOL,  0, 1, (void *)&DSIFarField,0, "retain only far-field contributions to DSIPFT"},
//
     {"MomentFile",PA_STRING,  1, 1, (void *)&MomentFile, 0, "name of dipole moment output file"},
     {"PSDFile",   PA_STRING,  1, 1, (void *)&PSDFile,    0, "name of panel source density file"},
     {"PlotSurfaceCurrents", PA_BOOL, 0, 1, (void *)&PlotSurfaceCurrents, 0,"generate surface current visualization files"},
     {"HDF5File",  PA_STRING,  1, 1, (void *)&HDF5File,   0, "name of HDF5 file for BEM matrix/vector export"},
//
     {"LogLevel",  PA_STRING,  1, 1, (void *)&LogLevel,   0, "none | terse | verbose | verbose2"},
     {"Cache",     PA_STRING,  1, 1, (void *)&Cache,      0, "read/write cache"},
     {"ReadCache", PA_STRING,  1, MAXCACHE,(void *)ReadCache,  &nReadCache, "read cache"},
     {"WriteCache",PA_STRING,  1, 1, (void *)&WriteCache, 0,                "write cache"},
//
     {0,0,0,0,0,0,0}
   };

  ProcessOptions(argc, argv, OSArray);
  if (GeoFile==0)
   OSUsage(argv[0],OSArray,"--geometry option is mandatory");

  /*******************************************************************/
  /* process frequency-related options to construct a list of        */
  /* frequencies at which to run calculations                        */
  /*******************************************************************/
  HVector *OmegaList1=0, *OmegaList2=0, *OmegaList=0;
  if (OmegaFile) // process --OmegaFile option if present
   { 
     OmegaList1=new HVector(OmegaFile,LHM_TEXT);
     if (OmegaList1->ErrMsg)
      ErrExit(OmegaList1->ErrMsg);
   }
  if (nOmegaVals>0) // process -- Omega options if present
   {
     OmegaList2=new HVector(nOmegaVals, LHM_COMPLEX);
     for(int n=0; n<nOmegaVals; n++)
      OmegaList2->SetEntry(n,OmegaVals[n]);
   }
  if (  OmegaList1 && !OmegaList2 )
   OmegaList=OmegaList1;
  else if ( !OmegaList1 && OmegaList2 )
   OmegaList=OmegaList2;
  else if (  OmegaList1 && OmegaList2  )
   OmegaList=Concat(OmegaList1, OmegaList2);
  else 
   OSUsage(argv[0], OSArray, "you must specify at least one frequency");

  /*******************************************************************/
  /* process incident-field-related options to construct the data    */
  /* used to generate the incident field in our scattering problem   */
  /*******************************************************************/
  if ( npwPol != npwDir )
   ErrExit("numbers of --pwPolarization and --pwDirection options must agree");
  if ( ngbPol != ngbDir || ngbDir!=ngbCenter || ngbCenter!=ngbWaist )
   ErrExit("numbers of --gbPolarization, --gbDirection, --gbCenter, and --gbWaist options must agree ");
  if ( npsLoc!=npsStrength )
   ErrExit("numbers of --psLocation and --psStrength options must agree");
  if ( nglbI0 != nglbWaist || nglbWaist!=nglbpolname || nglbpolname!=nglP || nglP!=nglL )
   ErrExit("numbers of --glbI0, --glbpolname, --P, --L, and --glbWaist options must agree ");

  IncField *IFDList=0, *IFD;
  int npw, ngb, nps, nglb, nghb; 
  for(npw=0; npw<npwPol; npw++)
   { IFD=new PlaneWave(pwPol + 3*npw, pwDir + 3*npw);
     IFD->Next=IFDList;
     IFDList=IFD;
   };
  for(ngb=0; ngb<ngbCenter; ngb++)
   { IFD=new GaussianBeam(gbCenter + 3*ngb, gbDir + 3*ngb, gbPol + 3*ngb, gbWaist[ngb]);
     IFD->Next=IFDList;
     IFDList=IFD;
   };
     GHBeam *GHBeamInit = 0; //new GHBeam(HM,HN); 
     GLBeam *GLBeamInit = 0; //new GLBeam(HM,HN);   
   for(nghb=0; nghb<ngHM; nghb++) 
    { printf("GHBeam being prepared\n");
      GHBeamInit = new GHBeam(HM[nghb],HN[nghb],ghbWaist[nghb]);
      IFD=GHBeamInit ; 
      IFD->Next=IFDList;
      IFDList=IFD;       
    }; 

   for(nglb=0; nglb<nglP; nglb++) 
    { printf("GLBeam being prepared\n");
      GLBeamInit = new GLBeam(P[nglb],L[nglb],glbWaist[nglb],glbI0[nglb],glbpolname[nglb]); 
      IFD=GLBeamInit ; 
      IFD->Next=IFDList;
      IFDList=IFD; 
    }; 
  for(nps=0; nps<npsLoc; nps++)
   { IFD=new PointSource(psLoc + 3*nps, psStrength + 3*nps);
     IFD->Next=IFDList;
     IFDList=IFD;
   };
  /*******************************************************************/
  /* sanity check to make sure the user specified an incident field  */
  /* if one is required for the outputs the user requested           */
  /*******************************************************************/
  bool NeedIncidentField = (    MomentFile!=0
                             || OPFTFile!=0
                             || EPPFTFile!=0
                             || DSIPFTFile!=0
                             || nEPFiles>0
                             || nFVMeshes>0
                             || PlotSurfaceCurrents
                           );
  if ( NeedIncidentField && IFDList==0 )
   ErrExit("you must specify at least one incident field source");

  //       printf("Sanity Checkpoint 1 past \n");
  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  SetLogFileName("OptTorque.log");
  Log("OptTorque running on %s",GetHostName());

  /*******************************************************************/
  /* PFT options *****************************************************/
  /*******************************************************************/
  PFTOptions MyPFTOpts, *PFTOpts=&MyPFTOpts;
  InitPFTOptions(PFTOpts);
  PFTOpts->DSIMesh     = DSIMesh;
  PFTOpts->DSIRadius   = DSIRadius;
  PFTOpts->DSIPoints   = DSIPoints;
  PFTOpts->DSIFarField = DSIFarField;
  PFTOpts->EPFTOrder   = EPFTOrder;
  //       printf("Sanity Checkpoint 2 past \n");
  /*******************************************************************/
  /* create the SSData structure containing everything we need to    */
  /* execute scattering calculations                                 */
  /*******************************************************************/
  SSData MySSData, *SSD=&MySSData;
  //       printf("Sanity Checkpoint 3 past \n");
  RWGGeometry *G = SSD->G = new RWGGeometry(GeoFile);
  HMatrix *M = SSD->M =SSD->G->AllocateBEMMatrix();
  SSD->RHS = SSD->G->AllocateRHSVector();
  SSD->IF  = IFDList;  
  HVector *KN = SSD->KN =SSD->G->AllocateRHSVector();
  //       printf("Sanity Checkpoint 4 past \n");
  char GeoFileBase[MAXSTR];
  strncpy(GeoFileBase, GetFileBase(GeoFile), MAXSTR);

  if (LogLevel) G->SetLogLevel(LogLevel);

  /*******************************************************************/
  /* sanity check: for now (20120924), calculations involving        */
  /* extended geometries must have only a single incident field      */
  /* source, which must be a plane wave, and the bloch wavevector    */
  /* is extracted from the plane wave direction                      */
  /*******************************************************************/
  double kBlochBuffer[3];
  if (G->LDim>0)
   { if ( npwPol!=1 || ngbCenter!=0 || npsLoc!=0 )
      ErrExit("for extended geometries, the incident field must be a single plane wave");
     SSD->kBloch = kBlochBuffer;
   }
  else
   SSD->kBloch=0;
  //       printf("Sanity Checkpoint 5 past \n");
  /*******************************************************************/
  /* preload the scuff cache with any cache preload files the user   */
  /* may have specified                                              */
  /*******************************************************************/
  if ( Cache!=0 && WriteCache!=0 )
   ErrExit("--cache and --writecache options are mutually exclusive");
  if (Cache) 
   WriteCache=Cache;
  for (int nrc=0; nrc<nReadCache; nrc++)
   PreloadCache( ReadCache[nrc] );
  if (Cache)
   PreloadCache( Cache );

  /*******************************************************************/
  /*******************************************************************/
  /*******************************************************************/
  void *HDF5Context=0;
  if (HDF5File)
   HDF5Context=HMatrix::OpenHDF5Context(HDF5File);

  //       printf("Sanity Checkpoint 6 past \n");
  //  double **ByEdge = (PlotPFTFlux ? AllocateByEdgeArray(SSD->G, 0) : 0);
  /*******************************************************************/
  char IFilename[MAXSTR];
  snprintf(IFilename,MAXSTR,"Intensity_%s.dat",GeoFileBase);  
  FILE *fIntensity=fopen(IFilename,"w"); 
  fprintf(fIntensity,"Omega, Intensity Integrated\n"); 
  fclose(fIntensity);

  /* loop over frequencies *******************************************/
  /*******************************************************************/
  char OmegaStr[MAXSTR];
  cdouble Omega;
  cdouble Eps, Mu;
  double wvnm;

  for(int nFreq=0; nFreq<OmegaList->N; nFreq++)
   { 
     //       printf("Sanity Checkpoint 7 past \n");
     Omega = OmegaList->GetEntry(nFreq);
     wvnm = 2.0*M_PI*1000.0/real(Omega); 
     z2s(Omega, OmegaStr);
     Log("Working at frequency %s...",OmegaStr);
     /*******************************************************************/
     /* assemble the BEM matrix at this frequency                       */
     /*******************************************************************/
     Log("Assembling BEM matrix...");
     if ( G->LDim==0 )
      G->AssembleBEMMatrix(Omega, M);
     else
      { cdouble EpsExterior, MuExterior;
        G->RegionMPs[0]->GetEpsMu(Omega, &EpsExterior, &MuExterior);
        double kExterior = real( csqrt2(EpsExterior*MuExterior) * Omega );
        SSD->kBloch[0] = kExterior*pwDir[0];
        SSD->kBloch[1] = kExterior*pwDir[1];
        SSD->kBloch[2] = 0.0;
        G->AssembleBEMMatrix(Omega, SSD->kBloch, M);
      }
     Log("Assembling RHS vector...");
     //         printf("Sanity Checkpoint 8 past \n");

     /*******************************************************************/
     /* dump the scuff cache to a cache storage file if requested. note */
     /* we do this only once per execution of the program, after the    */
     /* assembly of the BEM matrix at the first frequency, since at that*/
     /* point all cache elements that are to be computed will have been */
     /* computed and the cache will not grow any further for the rest   */
     /* of the program run.                                             */
     /*******************************************************************/
     if (WriteCache)
       { StoreCache( WriteCache );
         WriteCache=0;       
       }

     /*******************************************************************/
     /* export BEM matrix to a binary .hdf5 file if that was requested  */
     /*******************************************************************/
     if (HDF5Context)
       M->ExportToHDF5(HDF5Context,"M_%s",OmegaStr);

     /*******************************************************************/
     /* if the user requested no output options (for example, if she   **/
     /* just wanted to export the matrix to a binary file), don't      **/
     /* bother LU-factorizing the matrix or assembling the RHS vector. **/
     /*******************************************************************/
     if ( !NeedIncidentField )
       continue;

     /***************************************************************/
     /* set up the incident field profile and assemble the RHS vector */
     /***************************************************************/
     Log("  Assembling the RHS vector..."); 
     G->AssembleRHSVector(Omega, SSD->kBloch, IFDList, KN);

     double Intensity=GetIntegratedIntensity(G, 0, KN);
     printf("Integrated intensity=%e.\n",Intensity);
     fIntensity=fopen(IFilename,"a");
     fprintf(fIntensity,"%s   %e\n",OmegaStr,Intensity); 
     fclose(fIntensity); 
     SSD->RHS->Copy(SSD->KN); // copy RHS vector for later 

     /*******************************************************************/
     /* LU-factorize the BEM matrix to prepare for solving scattering   */
     /* problems                                                        */
     /*******************************************************************/
     Log("  LU-factorizing BEM matrix...");
     M->LUFactorize();
     //         printf("Sanity Checkpoint 9 past \n");
     /***************************************************************/
     /* solve the BEM system*****************************************/
     /***************************************************************/
     Log("  Solving the BEM system...");
     M->LUSolve(KN);

     if (HDF5Context)
      { SSD->RHS->ExportToHDF5(HDF5Context,"RHS_%s",OmegaStr);
        SSD->KN->ExportToHDF5(HDF5Context,"KN_%s",OmegaStr);
      }
     //         printf("Sanity Checkpoint 8 past \n");

     /*--------------------------------------------------------------*/
     /*- export the matrices into separate data files if asked  -----*/
     /*--------------------------------------------------------------*/
     if(MatrixOut)
       {
         char MatFileName[100];
         // store int(lambda[nm]) rather than omega in the filename. 
         snprintf(MatFileName, 100, "Mat_RHS_%inm.dat",int(wvnm));  
         SSD->RHS->ExportToText(MatFileName,"--separate,"); 
         snprintf(MatFileName, 100,"Mat_M_%inm.dat",int(wvnm));
         SSD->M->ExportToText(MatFileName,"--separate,"); 
         snprintf(MatFileName, 100,"Mat_KN_%inm.dat",int(wvnm));
         SSD->KN->ExportToText(MatFileName,"--separate,");  
       }

     /***************************************************************/
     /* now process all requested outputs                           */
     /***************************************************************/
     SSD->Omega=Omega;
     //         printf("Sanity Checkpoint 9 past \n");

     /*--------------------------------------------------------------*/
     /*- power, force, torque by various methods --------------------*/
     /*--------------------------------------------------------------*/
     if (OPFTFile)
      WritePFTFile(SSD, PFTOpts, SCUFF_PFT_OVERLAP, PlotPFTFlux, OPFTFile);
     //         printf("Sanity Checkpoint 10 past \n");
     if (EPPFTFile)
      WritePFTFile(SSD, PFTOpts, SCUFF_PFT_EP, PlotPFTFlux, EPPFTFile);
     //         printf("Sanity Checkpoint 11 past \n");
     if (DSIPFTFile)
      WritePFTFile(SSD, PFTOpts, SCUFF_PFT_DSI, PlotPFTFlux, DSIPFTFile);
     //         printf("Sanity Checkpoint 12 past \n");
     if (PFTFile) // default is overlap + EP for scattered power
      WritePFTFile(SSD, PFTOpts, SCUFF_PFT_DEFAULT, PlotPFTFlux, PFTFile);
     //         printf("Sanity Checkpoint 13 past \n");
     /*--------------------------------------------------------------*/
     /*- panel source densities -------------------------------------*/
     /*--------------------------------------------------------------*/
     if (PSDFile)
      WritePSDFile(SSD, PSDFile);
 
     /*--------------------------------------------------------------*/
     /*- scattered fields at user-specified points ------------------*/
     /*--------------------------------------------------------------*/
     int nepf;
     for(nepf=0; nepf<nEPFiles; nepf++)
      ProcessEPFile(SSD, EPFiles[nepf]);

     /*--------------------------------------------------------------*/
     /*- induced dipole moments       -------------------------------*/
     /*--------------------------------------------------------------*/
     if (MomentFile)
      GetMoments(SSD, MomentFile);

     /*--------------------------------------------------------------*/
     /*- surface current visualization-------------------------------*/
     /*--------------------------------------------------------------*/
     if (PlotSurfaceCurrents)
      G->PlotSurfaceCurrents(KN, Omega, "%s.%s.pp",GetFileBase(GeoFile),z2s(Omega));

     /*--------------------------------------------------------------*/
     /*- field visualization meshes ---------------------------------*/
     /*--------------------------------------------------------------*/
     int nfm;
     const char *FMMeshFilename; 
     for(nfm=0; nfm<nFVMeshes; nfm++){
       VisualizeFields(SSD, FVMeshes[nfm]);};
     //         printf("Sanity Checkpoint 12 past \n");
   }; //  for(nFreq=0; nFreq<NumFreqs; nFreqs++)
    
  /***************************************************************/
  /***************************************************************/
  /***************************************************************/
  if (HDF5Context)
   HMatrix::CloseHDF5Context(HDF5Context);

  //         printf("Sanity Checkpoint 13 past \n");
  delete M;
  delete KN;
  delete G;

  printf("Thank you for your support.\n");
}//end main 


void VisualizeFields(RWGGeometry *G, IncField *IF, HVector *KN,
                     cdouble Omega, char *MeshFileName)
{
  /*--------------------------------------------------------------*/
  /*- try to open output file ------------------------------------*/
  /*--------------------------------------------------------------*/
  char GeoFileBase[100], PPFileName[100];
  strncpy(GeoFileBase,GetFileBase(G->GeoFileName),100);
  snprintf(PPFileName,100,"%s.%s.pp",GeoFileBase,GetFileBase(MeshFileName));
  FILE *f=fopen(PPFileName,"a");
  if (!f)
   ErrExit("could not open field visualization file %s",PPFileName);

  /*--------------------------------------------------------------*/
  /*- try to open user's mesh file -------------------------------*/
  /*--------------------------------------------------------------*/
  RWGSurface *S=new RWGSurface(MeshFileName);

  Log("Creating flux plot for surface %s...",MeshFileName);
  printf("Creating flux plot for surface %s...\n",MeshFileName);

  /*--------------------------------------------------------------*/
  /*- create an Nx3 HMatrix whose columns are the coordinates of  */
  /*- the flux mesh panel vertices                                */
  /*--------------------------------------------------------------*/
  HMatrix *XMatrix=new HMatrix(S->NumVertices, 3);

  for(int nv=0; nv<S->NumVertices; nv++)
   {
     XMatrix->SetEntry(nv, 0, S->Vertices[3*nv + 0]);
     XMatrix->SetEntry(nv, 1, S->Vertices[3*nv + 1]);
     XMatrix->SetEntry(nv, 2, S->Vertices[3*nv + 2]);
   };

  /*--------------------------------------------------------------*/
  /*- get the total fields at the panel vertices                 -*/
  /*--------------------------------------------------------------*/
  HMatrix *FMatrix=G->GetFields(IF, KN, Omega, 0,
                                XMatrix, 0, FieldFuncs);

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  for(int nff=0; nff<NUMFIELDFUNCS; nff++)
   {
     fprintf(f,"View \"%s(%s)\" {\n",FieldTitles[nff],z2s(Omega));

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     for(int np=0; np<S->NumPanels; np++)
      {
        RWGPanel *P=S->Panels[np];
        int iV1 = P->VI[0];  double *V1 = S->Vertices + 3*iV1;
        int iV2 = P->VI[1];  double *V2 = S->Vertices + 3*iV2;
        int iV3 = P->VI[2];  double *V3 = S->Vertices + 3*iV3;

        fprintf(f,"ST(%e,%e,%e,%e,%e,%e,%e,%e,%e) {%e,%e,%e};\n",
                   V1[0], V1[1], V1[2],
                   V2[0], V2[1], V2[2],
                   V3[0], V3[1], V3[2],
                   FMatrix->GetEntryD(iV1,nff),
                   FMatrix->GetEntryD(iV2,nff),
                   FMatrix->GetEntryD(iV3,nff));
      };

     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     /*--------------------------------------------------------------*/
     fprintf(f,"};\n\n");
   };
  fclose(f);

  delete FMatrix;
  delete XMatrix;
  delete S;
}

/***************************************************************/
/***************************************************************/
/***************************************************************/
double **AllocateByEdgeArray(RWGGeometry *G, int ns)
{
  int NE = G->Surfaces[ns]->NumEdges;
  double **ByEdge=(double **)mallocEC(7*sizeof(double *));
  ByEdge[0]=(double *)mallocEC(7*NE*sizeof(double));
  for(int nq=1; nq<7; nq++)
   ByEdge[nq] = ByEdge[nq-1] + NE;

  return ByEdge;
}

void ProcessByEdgeArray(RWGGeometry *G, int ns, cdouble Omega,
                        double **ByEdge)
{
  static const char *PFTNames[7]
   ={"PAbs","FX","FY","FZ","TX","TY","TZ"};

  char FileName[100];
  snprintf(FileName,100,"%s.pp",GetFileBase(G->GeoFileName));

  for(int nq=0; nq<7; nq++)
   { char Tag[20];
     snprintf(Tag,20,"%s(%s)",PFTNames[nq],z2s(Omega));
     G->Surfaces[ns]->PlotScalarDensity(ByEdge[nq],true,FileName,Tag);
   };

 // free(ByEdge[0]);
 // free(ByEdge);

}
