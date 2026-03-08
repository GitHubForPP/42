/*    This file is distributed with 42,                               */
/*    the (mostly harmless) spacecraft dynamics simulation            */
/*    created by Eric Stoneking of NASA Goddard Space Flight Center   */

/*    Copyright 2010 United States Government                         */
/*    as represented by the Administrator                             */
/*    of the National Aeronautics and Space Administration.           */

/*    No copyright is claimed in the United States                    */
/*    under Title 17, U.S. Code.                                      */

/*    All Other Rights Reserved.                                      */


#include "42.h"

/* #ifdef __cplusplus
** namespace _42 {
** using namespace Kit;
** #endif
*/

/******************************************************************************/
/* References:                                                                */
/*                                                                            */
/* [1]  Louis J. Ippolito, Jr.                                                */
/*   "Satellite Communications Systems Engineering:                           */
/*   Atmospheric Effects, Satellite Link Design and System Performance"       */
/*   Wiley 2008.                                                              */
/*                                                                            */
/* [2]  Wertz et al. "Space Mission Engineering: The New SMAD",               */
/*   Microcosm 2011.                                                          */

/******************************************************************************/
long IsOcculted(double Pos[3], double Dir[3], double Range, double rad)
{
      double P,PosDir[3],PoD,RbyP;
      
      P = CopyUnitV(Pos,PosDir);
      PoD = VoV(PosDir,Dir);
      RbyP = rad/P;

      if (PoD < 0.0 && RbyP*RbyP > 1.0-PoD*PoD) return(1);
      else return(0);
}
/******************************************************************************/
void PosToENU(double PosW[3],double CPW[3][3])
{
      double Zaxis[3] = {0.0,0.0,1.0};
      double Up[3],North[3],East[3];
      long i;

      CopyUnitV(PosW,Up);
      VxV(Zaxis,Up,East);
      UNITV(East);
      VxV(Up,East,North);
      UNITV(North);
      for(i=0;i<3;i++) {
         CPW[0][i] = East[i];
         CPW[1][i] = North[i];
         CPW[2][i] = Up[i];
      }
}
/******************************************************************************/
void GndPosVel(struct WorldType *W, struct GroundStationType *G, double Delay, 
   double PosN[3], double VelN[3], double PosH[3], double VelH[3])
{
      double CDW[3][3]; /* Delay */
      double CDN[3][3];
      double Zaxis[3] = {0.0,0.0,1.0};
      long i;

      SimpRot(Zaxis,-Delay*W->w,CDW);
      MxM(CDW,W->CWN,CDN);
      MTxV(CDN,G->PosW,PosN);
      VelN[0] = -W->w*PosN[1];
      VelN[1] =  W->w*PosN[0];
      VelN[2] = 0.0; 

      MTxV(W->CNH,PosN,PosH);
      MTxV(W->CNH,VelN,VelH);
      for(i=0;i<3;i++) {
         PosH[i] += W->PosH[i];
         VelH[i] += W->VelH[i];
      }
}
/******************************************************************************/
void ScPosVel(struct WorldType *W, struct OrbitType *O, double Delay, 
   double PosN[3], double VelN[3], double PosH[3], double VelH[3])
{
      double anom;
      long i;
      
      //if (O->Regime == ORB_THREE_BODY) {
      //   if (O->LagDOF == LAGDOF_MODES) {
      //      LagModes2RV(Time,&LagSys[O->Sys],
      //         O,O->PosN,O->VelN);
      //   }
      //   else if (O->LagDOF == LAGDOF_COWELL) {
      //      ThreeBodyOrbitRK4(O);
      //      RV2LagModes(Time,&LagSys[O->Sys],O);
      //      O->Epoch = Time;
      //   }
      //   else if (O->LagDOF == LAGDOF_SPLINE) {
      //      SplineToPosVel(O);
      //   }
      //}
      //else if (O->Regime == ORB_CENTRAL) {
         //if (O->SplineActive) SplineToPosVel(O);
         //else if (O->J2DriftEnabled) MeanEph2RV(O,Time);
         //else {
            Eph2RV(O->mu,O->SLR,O->ecc,O->inc,O->RAAN,O->ArgP,
                   DynTime-O->tp-Delay,PosN,VelN,&anom);
         //}
      //}
      ///* ORB_FLIGHT is handled in Ephemerides() */
      ///* ORB_ZERO requires no action */

      MTxV(W->CNH,PosN,PosH);
      MTxV(W->CNH,VelN,VelH);
      for(i=0;i<3;i++) {
         PosH[i] += W->PosH[i];
         VelH[i] += W->VelH[i];
      }
}
/******************************************************************************/
/* Tx is Gnd, Rx is S/C                                                       */
void UplinkGeometry(struct CommLinkType *L)
{
      struct GroundStationType *G;
      struct SCType *S;
      struct WorldType *Wtx,*Wrx;
      struct OrbitType *O;
      struct BodyType *B;
      double CAW[3][3];

      G = &GroundStation[L->TxID];
      S = &SC[L->RxID];
      B = &S->B[L->RxBody];
      O = &Orb[S->RefOrb];
      L->TxWorld = G->World;
      L->RxWorld = O->World;
      Wtx = &World[L->TxWorld];
      Wrx = &World[L->RxWorld];

      GndPosVel(Wtx,G,L->Delay,L->TxPosN,L->TxVelN,L->TxPosH,L->TxVelH);
      /* Assume ground station "A" frame is ENU */
      PosToENU(G->PosW,CAW);
      MxM(CAW,Wtx->CWN,L->TxCAN);
      MxM(L->TxCAN,Wtx->CNH,L->TxCAH);

      ScPosVel(Wrx,O,L->Delay,L->RxPosN,L->RxVelN,L->RxPosH,L->RxVelH);
      MxM(L->RxCAB,B->CN,L->RxCAN);
      MxM(L->RxCAN,Wrx->CNH,L->RxCAH);
}
/******************************************************************************/
/* Tx is S/C, Rx is Gnd                                                       */
void DownlinkGeometry(struct CommLinkType *L)
{
      struct GroundStationType *G;
      struct SCType *S;
      struct WorldType *Wtx,*Wrx;
      struct OrbitType *O;
      struct BodyType *B;
      double CAW[3][3];

      G = &GroundStation[L->RxID];
      S = &SC[L->TxID];
      B = &S->B[L->TxBody];
      O = &Orb[S->RefOrb];
      L->TxWorld = O->World;
      L->RxWorld = G->World;
      Wtx = &World[L->TxWorld];
      Wrx = &World[L->RxWorld];

      GndPosVel(Wrx,G,L->Delay,L->RxPosN,L->RxVelN,L->RxPosH,L->RxVelH);
      /* Assume ground station "A" frame is ENU */
      PosToENU(G->PosW,CAW);
      MxM(CAW,Wrx->CWN,L->RxCAN);
      MxM(L->RxCAN,Wrx->CNH,L->RxCAH);

      ScPosVel(Wtx,O,L->Delay,L->TxPosN,L->TxVelN,L->TxPosH,L->TxVelH);
      MxM(L->TxCAB,B->CN,L->TxCAN);
      MxM(L->TxCAN,Wtx->CNH,L->TxCAH);
}
/******************************************************************************/
/* Both Tx and Rx are S/C                                                     */
void CrosslinkGeometry(struct CommLinkType *L)
{
      struct SCType *Stx,*Srx;
      struct WorldType *Wtx,*Wrx;
      struct OrbitType *Otx,*Orx;
      struct BodyType *Btx,*Brx;

      Stx = &SC[L->TxID];
      Btx = &Stx->B[L->TxBody];
      Otx = &Orb[Stx->RefOrb];
      Srx = &SC[L->RxID];
      Brx = &Srx->B[L->RxBody];
      Orx = &Orb[Srx->RefOrb];
      L->RxWorld = Orx->World;
      L->TxWorld = Otx->World;
      Wtx = &World[L->TxWorld];
      Wrx = &World[L->RxWorld];

      MxM(L->TxCAB,Btx->CN,L->TxCAN);
      MxM(L->TxCAN,Wtx->CNH,L->TxCAH);

      MxM(L->RxCAB,Brx->CN,L->RxCAN);
      MxM(L->RxCAN,Wrx->CNH,L->RxCAH);

      ScPosVel(Wtx,Otx,L->Delay,L->TxPosN,L->TxVelN,L->TxPosH,L->TxVelH);
      ScPosVel(Wrx,Orx,L->Delay,L->RxPosN,L->RxVelN,L->RxPosH,L->RxVelH);
}
/******************************************************************************/
void LinkGeometry(struct CommLinkType *L)
{
      struct WorldType *Wtx,*Wrx;
      double RelVelN[3],RelVelH[3];
      double dt;
      double PrevDelay = 0.0;
      long i;

      Wtx = &World[L->TxWorld];
      Wrx = &World[L->RxWorld];

      L->Delay = 0.0;
      do {
         if (L->LinkType == UPLINK) UplinkGeometry(L);
         else if (L->LinkType == DOWNLINK) DownlinkGeometry(L);
         else if (L->LinkType == CROSSLINK) CrosslinkGeometry(L);
         else {
            printf("Error: Unknown link type.\n");
            exit(1);
         }
   
         /* Find Range, Delay */
         if (L->TxWorld == L->RxWorld) {
            for(i=0;i<3;i++) {
               L->TxPathDirN[i] = L->RxPosN[i] - L->TxPosN[i];
            }
            L->Range = UNITV(L->TxPathDirN);
            L->Delay = L->Range/SPEED_OF_LIGHT;
         }
         else {
            for(i=0;i<3;i++) {
               L->TxPathDirH[i] = L->RxPosH[i] - L->TxPosH[i];
            }
            L->Range = UNITV(L->TxPathDirH);
            L->Delay = L->Range/SPEED_OF_LIGHT;
         }
         dt = L->Delay - PrevDelay;
         PrevDelay = L->Delay;
      } while(fabs(dt) > L->DelayTol);
      
      /* Find RangeRate, Doppler */
      if (L->TxWorld == L->RxWorld) {
         MTxV(Wtx->CNH,L->TxPathDirN,L->TxPathDirH);
         for(i=0;i<3;i++) {
            L->RxPathDirN[i] = -L->TxPathDirN[i];
            L->RxPathDirH[i] = -L->TxPathDirH[i];
            RelVelN[i] = L->TxVelN[i] - L->RxVelN[i];
         }
         L->RangeRate = VoV(RelVelN,L->RxPathDirN);
         L->Doppler = -L->RangeRate/SPEED_OF_LIGHT;
      }
      else {
         MxV(Wtx->CNH,L->TxPathDirH,L->TxPathDirN);
         for(i=0;i<3;i++) {
            L->RxPathDirH[i] = -L->TxPathDirH[i];
            RelVelH[i] = L->TxVelH[i] - L->RxVelH[i];
         }
         MxV(Wrx->CNH,L->RxPathDirH,L->RxPathDirN);
         L->RangeRate = VoV(RelVelH,L->RxPathDirH);
         L->Doppler = -L->RangeRate/SPEED_OF_LIGHT;
      }

      /* Check for occultation */
      if (L->LinkType == UPLINK) {
         if (VoV(L->TxPosN,L->TxPathDirN) > 0.0) L->TxOcculted = 0;
         else L->TxOcculted = 1;
         L->RxOcculted = IsOcculted(L->RxPosN,L->RxPathDirN,L->Range,Wrx->rad);
      }
      else if (L->LinkType == DOWNLINK) {
         L->TxOcculted = IsOcculted(L->TxPosN,L->TxPathDirN,L->Range,Wtx->rad);
         if (VoV(L->RxPosN,L->RxPathDirN) > 0.0) L->RxOcculted = 0;
         else L->RxOcculted = 1;
      }
      else {
         L->TxOcculted = IsOcculted(L->TxPosN,L->TxPathDirN,L->Range,Wtx->rad);
         L->RxOcculted = IsOcculted(L->RxPosN,L->RxPathDirN,L->Range,Wrx->rad);
      }
      if (L->TxOcculted || L->RxOcculted) L->PathIsOcculted = 1;
      else L->PathIsOcculted = 0;

}
/******************************************************************************/
/*  This is a simple stochastic model of losses through an atmosphere.        */
/*  Feel free to replace with your own loss model.                            */
void AtmoLossModel(struct CommLinkType *L)
{
      if (L->LinkType == CROSSLINK) L->AtmoPathLoss = 0.0;
      else {
         L->AtmoPathLoss = L->AtmoMean 
            + (1.0-DTSIM/L->AtmoCorrTime)*(L->AtmoPathLoss-L->AtmoMean) 
            + L->AtmoStd*GaussianRandom(L->AtmoNoise);
      }
}
/******************************************************************************/
void TxPerformance(struct CommLinkType *L)
{
      struct AntPattType *A;
      struct MeshType *M;
      double ViewPt[3],ViewDir[3];
      double TxPathDirA[3],ProjPt[3];
      long Ip,FoundPoly;
      long i;

      A = &L->TxAntPatt;
      M = &Mesh[A->MeshTag];

/* .. Antenna Gain */
      MxV(L->TxCAN,L->TxPathDirN,TxPathDirA);
      for(i=0;i<3;i++) {
         ViewPt[i] = 2.0*M->BBox.radius*TxPathDirA[i];
         ViewDir[i] = -TxPathDirA[i];
      }
      FoundPoly = OCProjectRayOntoMesh(ViewPt,ViewDir,
         M,ProjPt,&Ip);
      if (FoundPoly) {
         L->TxAntGain = (A->PeakGain-A->FloorGain)*MAGV(ProjPt)+A->FloorGain;
      }
      else {
         L->TxAntGain = A->FloorGain;
      }
}
/******************************************************************************/
void RxPerformance(struct CommLinkType *L)
{
      struct AntPattType *A;
      struct MeshType *M;
      double ViewPt[3],ViewDir[3];
      double RxPathDirA[3],ProjPt[3];
      long Ip,FoundPoly;
      long i;

      A = &L->RxAntPatt;
      M = &Mesh[A->MeshTag];

/* .. Antenna Gain */
      MxV(L->RxCAN,L->RxPathDirN,RxPathDirA);
      for(i=0;i<3;i++) {
         ViewPt[i] = 2.0*M->BBox.radius*RxPathDirA[i];
         ViewDir[i] = -RxPathDirA[i];
      }
      FoundPoly = OCProjectRayOntoMesh(ViewPt,ViewDir,
         M,ProjPt,&Ip);
      if (FoundPoly) {
         L->RxAntGain = (A->PeakGain-A->FloorGain)*MAGV(ProjPt)+A->FloorGain;
      }
      else {
         L->RxAntGain = A->FloorGain;
      }
}
/******************************************************************************/
void LoadLinkFile(void)
{
      FILE *infile;
      long Il,Seq;
      long EdgesEnabled = 1;
      double ang1,ang2,ang3;
      double RangeTol;
      double dt;
      char response[80],junk[80],newline;
      struct CommLinkType *L;
      long Seed;
      long OldNmesh;

      infile = FileOpen(InOutPath,"Inp_CommLink.txt","r");
      
      fscanf(infile,"%[^\n] %[\n]",junk,&newline);
      fscanf(infile,"%ld %[^\n] %[\n]",&Nlink,junk,&newline);
      
      CommLink = (struct CommLinkType *) calloc(Nlink,sizeof(struct CommLinkType));
      for(Il=0;Il<Nlink;Il++) {
         L = &CommLink[Il];
         fscanf(infile,"%[^\n] %[\n]",junk,&newline);  
         fscanf(infile,"%[^\n] %[\n]",junk,&newline); 

         fscanf(infile,"%s %lf  %[^\n] %[\n]",response,&dt,junk,&newline);
         L->OutEnabled = DecodeString(response);
         if (dt < DTSIM) {
            printf("Link[%ld] output timestep < DTSIM.  You'll want to fix that.\n",Il);
            exit(1);
         }
         else L->MaxOutCtr = (long) (dt/DTSIM+0.5);

         fscanf(infile,"%lf %[^\n] %[\n]",&RangeTol,junk,&newline); 
         L->DelayTol = RangeTol/SPEED_OF_LIGHT;

         fscanf(infile,"%s %[^\n] %[\n]",response,junk,&newline); 
         L->LinkType = DecodeString(response);
         fscanf(infile,"%lf %[^\n] %[\n]",&L->Freq,junk,&newline);
         L->Wavelength = SPEED_OF_LIGHT/L->Freq;
         fscanf(infile,"%ld %ld %[^\n] %[\n]",&L->TxID,&L->TxBody,junk,&newline);
         fscanf(infile,"%lf %lf %lf %ld %[^\n] %[\n]",
            &ang1,&ang2,&ang3,&Seq,junk,&newline);   
         A2C(Seq,ang1*D2R,ang2*D2R,ang3*D2R,L->TxCAB);      
         fscanf(infile,"%lf %[^\n] %[\n]",&L->TxPower,junk,&newline);
         fscanf(infile,"%lf %lf %[^\n] %[\n]",&L->TxAntPatt.PeakGain,
            &L->TxAntPatt.FloorGain,junk,&newline);
         fscanf(infile,"\"%[^\"]\" %[^\n] %[\n]",L->TxAntFileName,junk,&newline);
         fscanf(infile,"%ld %ld %[^\n] %[\n]",&L->RxID,&L->RxBody,junk,&newline);
         fscanf(infile,"%lf %lf %lf %ld %[^\n] %[\n]",
            &ang1,&ang2,&ang3,&Seq,junk,&newline);   
         A2C(Seq,ang1*D2R,ang2*D2R,ang3*D2R,L->RxCAB);
         fscanf(infile,"%lf %lf %[^\n] %[\n]",&L->RxAntPatt.PeakGain,
            &L->RxAntPatt.FloorGain,junk,&newline);
         fscanf(infile,"\"%[^\"]\" %[^\n] %[\n]",L->RxAntFileName,junk,&newline);
         fscanf(infile,"%lf %[^\n] %[\n]",&L->RxNoisePower,junk,&newline);
         fscanf(infile,"%lf %[^\n] %[\n]",&L->AtmoMean,junk,&newline);
         fscanf(infile,"%lf %lf %ld %[^\n] %[\n]",
            &L->AtmoStd,&L->AtmoCorrTime,&Seed,junk,&newline);
         L->AtmoNoise = CreateRandomProcess(Seed);
         
         OldNmesh = Nmesh;
         Mesh = LoadWingsObjFile(ModelPath,L->TxAntFileName,
            &Matl,&Nmatl,Mesh,&Nmesh,&L->TxAntPatt.MeshTag,EdgesEnabled);
         if (OldNmesh != Nmesh) LoadOctree(&Mesh[Nmesh-1]);
         OldNmesh = Nmesh;
         Mesh = LoadWingsObjFile(ModelPath,L->RxAntFileName,
            &Matl,&Nmatl,Mesh,&Nmesh,&L->RxAntPatt.MeshTag,EdgesEnabled);
         if (OldNmesh != Nmesh) LoadOctree(&Mesh[Nmesh-1]);
      }
      
      fclose(infile);
}
/******************************************************************************/
void CommLinkPerformance(void)
{
      struct CommLinkType *L;
      long Il;
      static long First = 1;

      if (First) {
         First = 0;
         LoadLinkFile();
      }

      for(Il=0;Il<Nlink;Il++) {
         L = &CommLink[Il];
         
         LinkGeometry(L);
         TxPerformance(L);
         RxPerformance(L);
         AtmoLossModel(L);
                  
         /* [1](4.6) */
         L->EIRP = L->TxPower + L->TxAntGain;
         
         /* [1](4.9) */
         L->PowerFluxDensity = L->EIRP - 20.0*log10(L->Range) - 10.99;
         
         /* [1](4.27) */
         L->FreeSpacePathLoss = 20.0*log10(TwoPi*L->Range/L->Wavelength);
         
         L->Loss = L->FreeSpacePathLoss + L->AtmoPathLoss;

         /* [1](4.31),(4.55) */
         L->Carrier = L->EIRP + L->RxAntGain - L->Loss;
         
         L->Noise = L->RxNoisePower;
                  
         L->CNR = L->Carrier - L->Noise;
      }
}

/* #ifdef __cplusplus
** }
** #endif
*/
