

/*
 * Include Files
 *
 */
#if defined(MATLAB_MEX_FILE)
#include "tmwtypes.h"
#include "simstruc_types.h"
#else
#include "rtwtypes.h"
#endif

/* %%%-SFUNWIZ_wrapper_includes_Changes_BEGIN --- EDIT HERE TO _END */
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "vsfmodel.h"
#define _USE_MATH_DEFINES
#include <math.h>
#include <simstruc.h>
#define FORCE_GAIN (1)
/* %%%-SFUNWIZ_wrapper_includes_Changes_END --- EDIT HERE TO _BEGIN */
#define u_width 1
#define y_width 1
/*
 * Create external references here.  
 *
 */
/* %%%-SFUNWIZ_wrapper_externs_Changes_BEGIN --- EDIT HERE TO _END */
/* extern double func(double a); */
/* %%%-SFUNWIZ_wrapper_externs_Changes_END --- EDIT HERE TO _BEGIN */

/*
 * Output functions
 *
 */
void vehmodelBUILDED_Outputs_wrapper(const real_T *inSteer,
			const real_T *inMz,
			real_T *outFz,
			real_T *outFy,
			real_T *outFx,
			real_T *outSlip,
			real_T *outStates,
			real_T *outDStates,
			const real_T *xD,
			const real_T *xC,
			const real_T  *params, const int_T  p_width0)
{
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_BEGIN --- EDIT HERE TO _END */
uint8_T i;


for( i=0; i<=4;i++)
{
    outStates[i] = xC[i];
    outDStates[i] = dstates[i];
}

for( i=FL;i<=RR;i++ )
{
    outFz[i] = xD[i+8];
    outFy[i] = xD[i+4];
    outFx[i] = xD[i];
    outSlip[i] = xD[i+12];
}
/* %%%-SFUNWIZ_wrapper_Outputs_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
  * Updates function
  *
  */
void vehmodelBUILDED_Update_wrapper(const real_T *inSteer,
			const real_T *inMz,
			const real_T *outFz,
			const real_T *outFy,
			const real_T *outFx,
			const real_T *outSlip,
			const real_T *outStates,
			const real_T *outDStates,
			real_T *xD,
			const real_T  *params,  const int_T  p_width0)
{
  /* %%%-SFUNWIZ_wrapper_Update_Changes_BEGIN --- EDIT HERE TO _END */
uint8_T i,j;
static bool isFirstIteration = 1;
real_T loadKN[4];
real_T slipDeg[4];
real_T mf_D,mf_E,mf_Bs;
real_T aux1,aux2,aux3;
const real_T aux_slip_angle_num[4] = {A,A,-B,-B};
const real_T aux_slip_angle_den[4] = {-TF2,TF2,-TR2,TR2};
const real_T aux_load_f1[4] = {ZF0,ZF0,ZR0,ZR0};
const real_T aux_load_f2[4] = {-Z[1],Z[1],-Z[6],Z[6]};
const real_T aux_load_f3[4] = {-Z[0],-Z[0],Z[0],Z[0]};
const real_T aux_load_f4[4] = {-Z[2],Z[2],-Z[4],Z[4]};
const real_T aux_load_f5[4] = {-Z[3],Z[3],-Z[5],Z[5]};

    if( isFirstIteration )
    {
        isFirstIteration = 0;
        ssPrintf("\n aceleracao inicial: %f , %f", LAT_ACCEL, LONG_ACCEL );
        for(j=0;j<=40;j++) 
        {
            ssPrintf("\n param %d %f", j, params[j]);
        }
    
    }

for( i=FL; i<= RR; i++ )
{

    slip_angle[i] = 0;
    if( ( i == FL) || ( i == FR) )
    {
        slip_angle[i] = inSteer[0]; 
    }
    
    slip_angle[i] -= atan2(outStates[0]+outStates[2]*aux_slip_angle_num[i],outStates[1]+outStates[2]*aux_slip_angle_den[i]);
    
    load[i] = aux_load_f1[i] +
              aux_load_f2[i]*LAT_ACCEL +
              aux_load_f3[i]*LONG_ACCEL +
              aux_load_f4[i]*outStates[4] +
              aux_load_f5[i]*outStates[3];

    loadKN[i] = load[i]*0.001;
    
    //aux1 = MFB[1]*loadKN[i] + MFB[2];
    //aux2 = ((MFB[3]*loadKN[i] + MFB[4])*exp(-MFB[5]*fabs(loadKN[i]))/(MFB[0]*aux1))*long_slip[i];
    //aux3 = loadKN[i]*(MFB[6]*loadKN[i] + MFB[7]) + MFB[8];
    //long_f[i] = loadKN[i]*aux1*sin(MFB[0]*atan((1-aux3)*aux2 + aux3*atan(aux2)));
    long_f[i] = 0;
    
    mf_D = (MFA[1]*loadKN[i] + MFA[2])*loadKN[i];
    mf_E = MFA[6]*loadKN[i] + MFA[7];
    mf_Bs = (MFA[3]*sin(2*atan(loadKN[i]/MFA[4]))/(MFA[0]*mf_D))*slip_angle[i]*180/M_PI;
    lat_f[i] = mf_D*sin(MFA[0]*atan((1-mf_E)*mf_Bs + mf_E*atan(mf_Bs)));      
    
    xD[i] = long_f[i];
    xD[i+4] = lat_f[i];
    xD[i+8] = load[i];
    xD[i+12]= slip_angle[i];
}
/* %%%-SFUNWIZ_wrapper_Update_Changes_END --- EDIT HERE TO _BEGIN */
}

/*
  *  Derivatives function
  *
  */
void vehmodelBUILDED_Derivatives_wrapper(const real_T *inSteer,
			const real_T *inMz,
			const real_T *outFz,
			const real_T *outFy,
			const real_T *outFx,
			const real_T *outSlip,
			const real_T *outStates,
			const real_T *outDStates,
			real_T *dx,
			real_T *xC,
			const real_T  *params,  const int_T  p_width0)
{
/* %%%-SFUNWIZ_wrapper_Derivatives_Changes_BEGIN --- EDIT HERE TO _END */
uint8_T i;
real_T aux_cos_steer,aux_sin_steer;
real_T veh_long_f[4], veh_lat_f[4];
real_T Mz,Fy,Fx;
const real_T aux_Mz_lat[4]= {A,A,-B,-B};
const real_T aux_Mz_long[4] = {-TF2,TF2,-TR2,TR2};

aux_cos_steer = cos(inSteer[0]);
aux_sin_steer = sin(inSteer[0]);
Mz = 0;
Fy = 0; 
Fx = 0;


for( i=FL; i<= RR; i++ )
{    
     
     if( ( i == FL) || ( i == FR) )
     {
         veh_long_f[i] = long_f[i]*aux_cos_steer - lat_f[i]*aux_sin_steer;
         veh_lat_f[i]  = long_f[i]*aux_sin_steer + lat_f[i]*aux_cos_steer;
     }
     else
     {
         veh_long_f[i] = long_f[i];
         veh_lat_f[i] = lat_f[i];                 
     }
    
    Mz+= aux_Mz_lat[i]*veh_lat_f[i] + aux_Mz_long[i]*veh_long_f[i];
    Fx += veh_long_f[i];
    Fy += veh_lat_f[i];
}

Mz += inMz[0];

dx[4] = xC[3];
dx[3] = D[0]*Fy + D[1]*Mz + D[2]*sin(xC[4]) - D[3]*xC[4] - D[4]*xC[3];
dx[2] = D[5]*Mz - D[6]*dx[3];
dx[1] = D[7]*Fx + D[8]*xC[2]*xC[0] - D[9]*xC[2]*xC[3];
dx[0] = D[10]*Fy + D[11]*dx[3] - D[12]*xC[2]*xC[1]; 

for( i=0; i<=4; i++)
    dstates[i] = dx[i];
/* %%%-SFUNWIZ_wrapper_Derivatives_Changes_END --- EDIT HERE TO _BEGIN */
}
