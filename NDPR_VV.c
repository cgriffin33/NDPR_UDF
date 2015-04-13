/**********************************************************************
   choud_pu.c                                                         
   UDF for Choudhuri Pitchup Case
   written by Christopher Griffin
   4-13-2015
***********************************************************************/

#include "udf.h"

/* Define constants in SI units */
#define RE 10000                                      /* Reynolds number */
#define VISC 1.7894E-5                                /* viscosity, N*s/m^2 */ 
#define DENS 1.225                                    /* density, kg/m^3 */
#define CHORD 0.2921                                  /* chord length, m */
#define NDPR 0.2                                      /* Non-dimensional pitch rate */

/* Define global variables */
static real aoa, aoa_old, prate, aprate, Vmag, t_o, t;

int open()
{
   /* open debug file for writing */
   FILE * fp;
   fp = fopen ("aoahistory.txt", "w+");
   fprintf(fp, "%s %s %s", "Time", "Old AOA", "Current AOA");
   fclose(fp);
   
   return (0);
}

DEFINE_PROFILE(x_velocity, thread, position) 
{
   face_t f;
   
   Vmag = (RE*VISC)/(DENS*CHORD); /*Velocity magnitude*/
   aprate = (NDPR*Vmag)/CHORD; /*Asymptotic pitching rate*/
   t_o = (0.5*CHORD)/Vmag; /*Time at which the pitch rate has reach 99% of asymptotic pitching rate*/
   t = CURRENT_TIME;
   
   /* Initialize aoa_old if this is first time step. */
   if (t == 0)
      {aoa_old = 0;}
   
   /* Define the variable pitching rate in rad/s */
   prate = aprate*(1-exp((-4.6*t)/t_o));
   
   /* Calculate current aoa */
   aoa = aoa_old + (t*prate);
   
   /* Loop through the inlet boundary and assign x velocity component */
   begin_f_loop(f, thread)
      {
         F_PROFILE(f, thread, position) = sqrt((Vmag*Vmag)/(1+(tan(aoa))));
      }
   end_f_loop(f, thread)
}

DEFINE_PROFILE(y_velocity, thread, position)
{
   face_t f;
   
   /* Loop through the inlet boundary and assign y velocity component */
   begin_f_loop(f, thread)
      {   
         F_PROFILE(f, thread, position) = tan(aoa)*sqrt((Vmag*Vmag)/(1+(tan(aoa))));
      }
   end_f_loop(f, thread)
   
   FILE * fp;
   fp = fopen ("aoahistory.txt", "w+");
   fprintf(fp, "%f.8 %f.8 %f.8", t, aoa_old, aoa);
   fclose(fp);
   
   /* current aoa becomes aoa_old */
   aoa_old = aoa;
}

