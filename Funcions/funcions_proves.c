#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "nr.h"
#include "nrutil.h"
#include "nrutil.c"
//#include "vectors.h"
//#include "vectors.c"
#include "tridag.c"
#include "cyclic.c"





/* ######################################################################### */

void  steric_force( vec_s  r0f, 
                    float  a,
                    vec_s *f_steric)
/* ######################################################################### */
{
    float d2, inv2, inv4, inv8, inv14;
    vec_s res;
    d2 = dotProduct(r0f,r0f);
    inv2 = a*a/d2;
    inv4 = inv2*inv2;
    inv8 = inv4*inv4;
    inv14 = inv8*inv4*inv2;
    res = scalarProduct(r0f, inv14);

    if(inv2>1.0){
      f_steric[0].x += res.x; 
      f_steric[0].y += res.y;
      f_steric[0].z += res.z;
    }
    else{
      f_steric[0].x += 0.0; 
      f_steric[0].y += 0.0;
      f_steric[0].z += 0.0;      
    }
    //printf("dins steric %le\t%le\t%le\t%le\t%le\t%le\n", d2, a, inv8, f_steric[0].x, f_steric[0].y,f_steric[0].z);

 
  }

/* ######################################################################### */

void  steric_force2( vec_s  r0f, 
                    float  a,
                    vec_s *f_steric)
/* ######################################################################### */
{
    float d2, inv2, sc1, inv6, inv12;
    a = a/pow(2.,1./6.);
    vec_s res;
    d2 = dotProduct(r0f,r0f);
    inv2 = a*a/d2;
    inv6 = inv2*inv2*inv2;
    inv12 = inv6*inv6;
    sc1 = (12.*inv12 - 6.*inv6)/d2;
    res = scalarProduct(r0f, sc1);

    if(inv6>0.5){
      f_steric[0].x += res.x; 
      f_steric[0].y += res.y;
      f_steric[0].z += res.z;
    }
    else{
      f_steric[0].x += 0.0; 
      f_steric[0].y += 0.0;
      f_steric[0].z += 0.0;      
    }
    //printf("dins steric %le\t%le\t%le\t%le\t%le\t%le\n", d2, a, inv8, f_steric[0].x, f_steric[0].y,f_steric[0].z);

 
  }


/* ######################################################################### */

float  bending_forces_1filament( vec_s  *rij1,
           float  *mags1,
           vec_s  *f,
           param_s params    )
/* ######################################################################### */
{


    int    n_mon1, i;
   
    float  b_l, k, phi_0, u, m_inv;
 
    vec_s  f_di[3];
  

    n_mon1 = params.n_mon1;
    b_l   = params.b_l;
    k     = params.k1;
    phi_0 = params.phi_0;
    u     = 0.0;
    m_inv = 1.0/params.mass1;
    for( i = 0; i < n_mon1 - 2; i++)
    {
       u += angle_force( rij1[i + 1], 
                         rij1[i+2], 
                         mags1[i+1], 
                         mags1[i+2],
                         k,
                         phi_0,
                         f_di);

       f[i].x += f_di[0].x*m_inv;
       f[i].y += f_di[0].y*m_inv;
       f[i].z += f_di[0].z*m_inv;
       f[i+1].x += f_di[1].x*m_inv;
       f[i+1].y += f_di[1].y*m_inv;
       f[i+1].z += f_di[1].z*m_inv;
       f[i+2].x += f_di[2].x*m_inv;
       f[i+2].y += f_di[2].y*m_inv;
       f[i+2].z += f_di[2].z*m_inv;
    }
    return( u );  
}

/* ########################################################################## */
float   angle_force( vec_s   r1,
                    vec_s   r2, 
                    float   r12,
                    float   r23,
                    float   k,
                    float   phi_0,
                    vec_s  *f )
/* ########################################################################## */
{
   int   i;

   float dprod, fac, u;
   float r12_3, r23_3, cphi, sphi, arg;

   dprod = r1.x*r2.x + r1.y*r2.y + r1.z*r2.z;
   cphi   = dprod/(r12*r23);

   u      = k*(1.0 - cphi);
   fac    = k;

   /*     sphi   = (1.0 - sqrt(cphi))/2.0; */
   /*     u      = k*(sphi)/2.0; */
   /*     fac    = k/4.0; */
 
   r12_3  = r12*r12*r12;
   r23_3  = r23*r23*r23;


   f[0].x =  r2.x/(r12*r23) - r1.x*dprod/(r12_3*r23);
   f[0].y =  r2.y/(r12*r23) - r1.y*dprod/(r12_3*r23);
   f[0].z =  r2.z/(r12*r23) - r1.z*dprod/(r12_3*r23);
   f[1].x = -r2.x*dprod/(r12*r23_3) + 
            (r1.x - r2.x)/(r12*r12) +
             r1.x*dprod/(r12_3*r23);
   f[1].y = -r2.y*dprod/(r12*r23_3) + 
            (r1.y - r2.y)/(r12*r12) +
             r1.y*dprod/(r12_3*r23);
   f[1].z = -r2.z*dprod/(r12*r23_3) + 
            (r1.z - r2.z)/(r12*r12) +
             r1.z*dprod/(r12_3*r23);
   f[2].x  = r2.x*dprod/(r12*r23_3) - r1.x/(r12*r23); 
   f[2].y  = r2.y*dprod/(r12*r23_3) - r1.y/(r12*r23); 
   f[2].z  = r2.z*dprod/(r12*r23_3) - r1.z/(r12*r23);
   for( i = 0; i < 3; i++ )
   {
      f[i].x *= fac;
      f[i].y *= fac;
      f[i].z *= fac;
   }
   return( u );
  
}


/* ########################################################################## */
void     constrain_positions_1filament( vec_s  *r1,
                              vec_s  *rij1,
                              vec_s  *v1,
                              vec_s  *fc_2,
                              param_s params )
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans;

    int   n_mon, n_c, i;
    float rijx, rijy, rijz, mag;
    float rhs_old, b_l, b_l_sq, max_error, error, tol;
    float dx, dy, dz, dt;

    n_mon   = params.n_mon1;
    n_c     = n_mon - 1;
    tol     = params.tolerance;
    dt      = params.dt;
    b_l     = params.b_l;
    b_l_sq  = b_l*b_l;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_2[i].x = 0.0;
        fc_2[i].y = 0.0;
        fc_2[i].z = 0.0;
    }
    ans[0]     = 0.0;
    ans[n_mon] = 0.0;
    rij1[0].x   = 0.0;
    rij1[0].y   = 0.0;
    rij1[0].z   = 0.0;
    rij1[n_c + 1].x = 0.0;
    rij1[n_c + 1].y = 0.0;
    rij1[n_c + 1].z = 0.0;

    rijx  = r1[0].x - r1[1].x;
    rijy  = r1[0].y - r1[1].y;
    rijz  = r1[0].z - r1[1].z;
    mag   = rijx*rijx + rijy*rijy + rijz*rijz; 

    a[1]   = 0.0;
    b[1]   = rij1[1].x*rijx;
    b[1]  += rij1[1].y*rijy;
    b[1]  += rij1[1].z*rijz;
    b[1]  *= 4.0;
    c[1]   = rij1[2].x*rijx;
    c[1]  += rij1[2].y*rijy;
    c[1]  += rij1[2].z*rijz;
    c[1]  *= -2.0;
    rhs[1] = b_l_sq - mag;

    rijx    = r1[n_c - 1].x - r1[n_c].x;
    rijy    = r1[n_c - 1].y - r1[n_c].y;
    rijz    = r1[n_c - 1].z - r1[n_c].z;
    mag     = rijx*rijx + rijy*rijy + rijz*rijz; 
    a[n_c]  = rij1[n_c - 1].x*rijx;
    a[n_c] += rij1[n_c - 1].y*rijy;
    a[n_c] += rij1[n_c - 1].z*rijz;
    a[n_c] *= -2.0;
    b[n_c]  = rij1[n_c].x*rijx;
    b[n_c] += rij1[n_c].y*rijy;
    b[n_c] += rij1[n_c].z*rijz;
    b[n_c] *= 4.0;
    c[n_c]  = 0.0;
    rhs[n_c]= b_l_sq - mag;
    for( i = 2; i < n_c; i++ )
    {
        rijx   = r1[i-1].x - r1[i].x;
        rijy   = r1[i-1].y - r1[i].y;
        rijz   = r1[i-1].z - r1[i].z;
        mag    = rijx*rijx + rijy*rijy + rijz*rijz; 
        a[i]   = -2.0*(rijx*rij1[i-1].x + rijy*rij1[i-1].y + rijz*rij1[i-1].z);
        b[i]   =  4.0*(rijx*rij1[i].x + rijy*rij1[i].y + rijz*rij1[i].z);
        c[i]   = -2.0*(rijx*rij1[i+1].x + rijy*rij1[i+1].y + rijz*rij1[i+1].z);
        rhs[i] = b_l_sq - mag;
    }
           //printf(" a constrain_positions_2filaments, abans de tridag \n");

    do
    { 
           //printf(" a constrain_positions_2filaments, abans de tridag \n");

       tridag( a, b, c, rhs, ans, n_c );
       //printf(" a constrain_positions_2filaments, despres de tridag \n");
       max_error = 0.0;
       for( i = 1; i <= n_c; i++ )
       {
          rijx    = r1[i-1].x + rij1[i].x*ans[i] - rij1[i-1].x*ans[i-1] -
                   (r1[i].x + rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i]);
          rijy    = r1[i-1].y + rij1[i].y*ans[i] - rij1[i-1].y*ans[i-1] -
                   (r1[i].y + rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i]);
          rijz    = r1[i-1].z + rij1[i].z*ans[i] - rij1[i-1].z*ans[i-1] -
                   (r1[i].z + rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i]);
          mag     = rijx*rijx + rijy*rijy + rijz*rijz;
          error   = (sqrt(mag) - b_l)/b_l;
          rhs_old = rhs[i];
          rhs[i]  = a[i]*ans[i-1] + b[i]*ans[i] + c[i]*ans[i+1];
          rhs[i]  = b_l_sq - mag + rhs[i];
          if( error > max_error )
          {
              max_error = error;
          }
       }
    }while( max_error > tol ); 

               //printf(" a constrain_positions_2filaments, despres de tridag \n");

    for( i = 0; i < n_mon; i++ )
    {
       dx      = rij1[i+1].x*ans[i+1] - rij1[i].x*ans[i];
       dy      = rij1[i+1].y*ans[i+1] - rij1[i].y*ans[i];
       dz      = rij1[i+1].z*ans[i+1] - rij1[i].z*ans[i];

       r1[i].x += dx;
       r1[i].y += dy;
       r1[i].z += dz;
     

       fc_2[i].x = 2.0*dx/(dt*dt);
       fc_2[i].y = 2.0*dy/(dt*dt);
       fc_2[i].z = 2.0*dz/(dt*dt);
    }
 
}




/* ########################################################################## */
void     constrain_velocities_1filament(  vec_s  *rij,
                                float  *mag,
                                vec_s  *v1,
        vec_s  *fc_1,
                                param_s params )
       
/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;

    n_mon   = params.n_mon1;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        b      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        c      = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float )); 
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_1[i].x = 0.0;
        fc_1[i].y = 0.0;
        fc_1[i].z = 0.0;
    }
    for( i = 1; i <= n_c; i++ )
    {
        rhs[i]  = (v1[i-1].x - v1[i].x)*rij[i].x;
        rhs[i] += (v1[i-1].y - v1[i].y)*rij[i].y;
        rhs[i] += (v1[i-1].z - v1[i].z)*rij[i].z;
        rhs[i] /= mag[i];
        a[i]    = rij[i-1].x*rij[i].x +
                  rij[i-1].y*rij[i].y +
                  rij[i-1].z*rij[i].z;
        b[i]    = -2.0*( rij[i].x*rij[i].x +
                   rij[i].y*rij[i].y +
                   rij[i].z*rij[i].z);
        c[i]    = rij[i].x*rij[i+1].x +
                  rij[i].y*rij[i+1].y +
                  rij[i].z*rij[i+1].z;
         a[i]  /= (mag[i]*mag[i]);
         b[i]  /= (mag[i]*mag[i]);
         c[i]  /= (mag[i]*mag[i]);
     }
     tridag( a, b, c, rhs, ans, n_c );
     for( i = 1; i <= n_c; i++ )
     {

       v1[i-1].x += ans[i]*rij[i].x/mag[i];
       v1[i-1].y += ans[i]*rij[i].y/mag[i];
       v1[i-1].z += ans[i]*rij[i].z/mag[i];
       v1[i].x   -= ans[i]*rij[i].x/mag[i];
       v1[i].y   -= ans[i]*rij[i].y/mag[i];
       v1[i].z   -= ans[i]*rij[i].z/mag[i];
       fc_1[i-1].x += ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i-1].y += ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i-1].z += ans[i]*rij[i].z/(dt*mag[i]);
       fc_1[i].x   -= ans[i]*rij[i].x/(dt*mag[i]);
       fc_1[i].y   -= ans[i]*rij[i].y/(dt*mag[i]);
       fc_1[i].z   -= ans[i]*rij[i].z/(dt*mag[i]);
       
     }
     


}

/* ########################################################################## */
void     constrain_forces_1filament(      vec_s  *rij1,
                                float  *mags1,
                                vec_s  *f,
                                vec_s  *fc_3,
                                param_s params )

/* ########################################################################## */
{
static int    init_flag = TRUE;
static float *a, *b, *c, *rhs, *ans, dt;

    int   n_mon, n_c, i;

    n_mon   = params.n_mon1;
    n_c     = n_mon - 1;

    if( init_flag )
    {
        a      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        b      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        c      = (float *)malloc( (n_mon + 2)*sizeof( float ));
        ans    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        rhs    = (float *)malloc( (n_mon + 2)*sizeof( float ));
        dt     = params.dt;
        init_flag = FALSE;
    }
    for( i = 0; i < n_mon; i++ )
    {
        fc_3[i].x = 0.0;
        fc_3[i].y = 0.0;
        fc_3[i].z = 0.0;
    }
    for( i = 1; i <= n_c; i++ )
    {
        rhs[i]  = (f[i-1].x - f[i].x)*rij1[i].x;
        rhs[i] += (f[i-1].y - f[i].y)*rij1[i].y;
        rhs[i] += (f[i-1].z - f[i].z)*rij1[i].z;
        rhs[i] /= mags1[i];
        a[i]    = rij1[i-1].x*rij1[i].x +
                  rij1[i-1].y*rij1[i].y +
                  rij1[i-1].z*rij1[i].z;
        b[i]    = -2.0*( rij1[i].x*rij1[i].x +
                   rij1[i].y*rij1[i].y +
                   rij1[i].z*rij1[i].z);
        c[i]    = rij1[i].x*rij1[i+1].x +
                  rij1[i].y*rij1[i+1].y +
                  rij1[i].z*rij1[i+1].z;
         a[i]  /= (mags1[i]*mags1[i]);
         b[i]  /= (mags1[i]*mags1[i]);
         c[i]  /= (mags1[i]*mags1[i]);
     }
     tridag( a, b, c, rhs, ans, n_c );
     for( i = 1; i <= n_c; i++ )
     {
         fc_3[i-1].x += ans[i]*rij1[i].x/(mags1[i]);
         fc_3[i-1].y += ans[i]*rij1[i].y/(mags1[i]);
         fc_3[i-1].z += ans[i]*rij1[i].z/(mags1[i]);
         fc_3[i].x   -= ans[i]*rij1[i].x/(mags1[i]);
         fc_3[i].y   -= ans[i]*rij1[i].y/(mags1[i]);
         fc_3[i].z   -= ans[i]*rij1[i].z/(mags1[i]);

     }
}


/*  ************************************************************************ */

void     oseen_1filament( vec_s *r,
    vec_s *fa,
    vec_s *v_h,
    param_s parameters)
/*  ************************************************************************ */
{
  static float hydr_radius, eta, gamma, bond, massa, massa1, massa2;
  static char  init_flag = TRUE;
  static int n_mon1, n_mon2, n_fil;
  
  float rijx, rijy,rijz;
  float mag_ij;
  int  i,j,k; 
  vec_s r_i, r_j, f_i, f_j, v_hydr;

  if( init_flag )
    {
      eta            = parameters.viscosity;
      gamma          = parameters.gamma1*parameters.mass1;
      massa          = parameters.mass1;
      n_mon1          = parameters.n_mon1;        
      init_flag  = FALSE;
    }

  
  for( i = 0; i < n_mon1; i++ )
    {
      v_hydr.x = 0;
      v_hydr.y = 0;
      v_hydr.z = 0;
      
      r_i.x      = r[i].x;
      r_i.y      = r[i].y;
      r_i.z      = r[i].z; 
      

      for( j = 0; j < n_mon1; j++)
  {
    

    f_j.x      = massa * fa[j].x;
    f_j.y      = massa * fa[j].y;
    f_j.z      = massa * fa[j].z;
    
    r_j.x      = r[j].x;
    r_j.y      = r[j].y;
    r_j.z      = r[j].z;
  
    rijx     = r[i].x - r[j].x;
    rijy     = r[i].y - r[j].y;
    rijz     = r[i].z - r[j].z;

    mag_ij  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
    
    rijx    /= mag_ij;
    rijy    /= mag_ij;
    rijz    /= mag_ij;
    
    
    if(j != i)
      {
        v_hydr.x += ((1+ rijx * rijx)* f_j.x +
         (rijx * rijy)   * f_j.y +
         (rijx * rijz)   * f_j.z  ) / mag_ij;
        
        v_hydr.y += ((rijx* rijy)    * f_j.x +
         (1+ rijy * rijy)* f_j.y +
         (rijy * rijz)   * f_j.z  ) / mag_ij;
    
        v_hydr.z += ((rijz * rijx)   * f_j.x +
         (rijz * rijy)   * f_j.y +
         (1+ rijz * rijz)* f_j.z  ) / mag_ij;
        if(i>n_mon1-2){
                }
      } 
  

  }

      v_h[i].x = ( v_hydr.x ) *  (1./(8.*PI*eta));
      v_h[i].y = ( v_hydr.y ) *  (1./(8.*PI*eta));
      v_h[i].z = ( v_hydr.z ) *  (1./(8.*PI*eta));


    }
}


/*  ************************************************************************ */

void     blake_1filament( vec_s *r,
    vec_s *fa,
    vec_s *v_h,
    param_s parameters)
/*  ************************************************************************ */
{
  static float hydr_radius, eta, gamma, bond, massa, massa1, massa2;
  static char  init_flag = TRUE;
  static int n_mon1, n_mon2, n_fil;
  
  float rijx, rijy,rijz,rijxbar, rijybar,rijzbar;
  float mag_ij, mag_ijbar, mag_ijbar2, mag_ijbar3, mag_ijbar5;
  int  i,j,k; 
  vec_s r_i, r_j, f_i, f_j, v_hydr;

  if( init_flag )
    {
      eta            = parameters.viscosity;
      gamma          = parameters.gamma1*parameters.mass1;
      massa          = parameters.mass1;
      n_mon1          = parameters.n_mon1;        
      init_flag  = FALSE;
    }



  /*    //printf("hydrodynamic radius is %f", hydr_radius); */
  
  for( i = 0; i < n_mon1; i++ )
    {

      v_hydr.x = 0;
      v_hydr.y = 0;
      v_hydr.z = 0;
      
      
      r_i.x      = r[i].x;
      r_i.y      = r[i].y;
      r_i.z      = r[i].z; 




      for( j = 0; j < n_mon1; j++)
  {
    

    f_j.x      = massa * fa[j].x;
    f_j.y      = massa * fa[j].y;
    f_j.z      = massa * fa[j].z;
    
    r_j.x      = r[j].x;
    r_j.y      = r[j].y;
    r_j.z      = r[j].z;
  
    rijx     = r[i].x - r[j].x;
    rijy     = r[i].y - r[j].y;
    rijz     = r[i].z - r[j].z;
    rijxbar     = r[i].x - r[j].x;
    rijybar     = r[i].y - r[j].y;
    rijzbar     = r[i].z + r[j].z;

    mag_ij  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
    mag_ijbar2  = rijx*rijx + rijy*rijy + rijzbar*rijzbar;
    mag_ijbar  = sqrt(mag_ijbar2);
    mag_ijbar3 = mag_ijbar*mag_ijbar2;
    mag_ijbar5 = mag_ijbar3*mag_ijbar2;
    
    rijx    /= mag_ij;
    rijy    /= mag_ij;
    rijz    /= mag_ij;
    rijxbar  /= mag_ijbar;
    rijybar  /= mag_ijbar;
    rijzbar  /= mag_ijbar;
    
    
    
    if(j != i)
      {
      // OSEEN:
        v_hydr.x +=  ((1+ rijx * rijx)* f_j.x +
         (rijx * rijy)   * f_j.y +
         (rijx * rijz)   * f_j.z  ) / mag_ij;
        
        v_hydr.y +=  ((rijx* rijy)    * f_j.x +
         (1+ rijy * rijy)* f_j.y +
         (rijy * rijz)   * f_j.z  ) / mag_ij;
    
        v_hydr.z +=  ((rijz * rijx)   * f_j.x +
         (rijz * rijy)   * f_j.y +
         (1+ rijz * rijz)* f_j.z  ) / mag_ij;
         
      // OSEEN_BAR:
        v_hydr.x -= ((1+ rijxbar * rijxbar)* f_j.x +
         (rijxbar * rijybar)   * f_j.y +
         (rijxbar * rijzbar)   * f_j.z  ) / mag_ijbar;
        
        v_hydr.y -= ((rijxbar* rijybar)    * f_j.x +
         (1+ rijybar * rijybar)* f_j.y +
         (rijybar * rijzbar)   * f_j.z  ) / mag_ijbar;
    
        v_hydr.z -= ((rijzbar * rijxbar)   * f_j.x +
         (rijzbar * rijybar)   * f_j.y +
         (1+ rijzbar * rijzbar)* f_j.z  ) / mag_ijbar;

      // deltaGim:
        
        v_hydr.x += -2*r_i.z*r_j.z*(1/mag_ijbar3 - 3*(r_i.x - r_j.x)*(r_i.x - r_j.x)/mag_ijbar5)* f_j.x + 
        6*(r_i.z*r_j.z*(r_i.x - r_j.x)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.y + 
        2*(r_i.x - r_j.x)*(r_j.z/mag_ijbar3 - 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
        
        v_hydr.y += 6*(r_i.z*r_j.z*(r_i.x - r_j.x)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.x -
        2*r_i.z*r_j.z*(1/mag_ijbar3 - 3*(r_i.y - r_j.y)*(r_i.y - r_j.y)/mag_ijbar5)* f_j.y + 
        2*(r_i.y - r_j.y)*(r_j.z/mag_ijbar3 - 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;
        
    
        v_hydr.z += 2*(r_i.x - r_j.x)*(r_j.z/mag_ijbar3 + 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.x +
        2*(r_i.y - r_j.y)*(r_j.z/mag_ijbar3 + 3*(r_i.z*r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.y +
        2*r_i.z*r_j.z*(1/mag_ijbar3- 3*(r_i.z + r_j.z)*(r_i.z + r_j.z)/mag_ijbar5)* f_j.z;

      }
    
    
  }

      v_h[i].x = ( v_hydr.x ) *  (1./(8.*PI*eta));
      v_h[i].y = ( v_hydr.y ) *  (1./(8.*PI*eta));
      v_h[i].z = ( v_hydr.z ) *  (1./(8.*PI*eta));
  
    }




}


/* ########################################################################## */
float    verlet_pt1_1filament(  vec_s  *r,
                      vec_s  *v,
                      vec_s  *fa,
                      vec_s  *vhydro,
                      param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq;
  static float hydr_radius, gamma;
  static float disp1;
  static float bond;
  static vec_s *v_h = NULL;
  static int  n_mon1;

  float ddx, ddy, ddz, f_filx;
  int  i;
  vec_s f_i, v_norm, v_hydr;
  
    if( init_flag )
    {

        dt             = parameters.dt;
  bond           = parameters.b_l;
  hydr_radius    = parameters.hydr_radius;
  n_mon1          = parameters.n_mon1;
  gamma = parameters.gamma1;


        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
        init_flag  = FALSE;
    }

    if(v_h == NULL){
      v_h = (vec_s *)malloc(n_mon1*sizeof( vec_s ));     
    }

    if (parameters.oseen ==1){
      oseen_1filament(r, fa, v_h, parameters);
    }
    if (parameters.blake ==1){
      blake_1filament(r, fa, v_h, parameters);
    }

    for( i = 0; i < n_mon1; i++ )
    {

        disp1         = 1.0/(1.0 + 0.5*gamma*dt);

        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;

  
        v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
    
        v_hydr.x   = v_h[i].x; 
        v_hydr.y   = v_h[i].y;
        v_hydr.z   = v_h[i].z;

        vhydro[i].x   = v_h[i].x; 
        vhydro[i].y   = v_h[i].y;
        vhydro[i].z   = v_h[i].z;

        ddx     = v_norm.x*dt + (f_i.x + gamma*(v_hydr.x-v_norm.x))*half_dt_sq;
        ddy     = v_norm.y*dt + (f_i.y + gamma*(v_hydr.y-v_norm.y))*half_dt_sq;
        ddz     = v_norm.z*dt + (f_i.z + gamma*(v_hydr.z-v_norm.z))*half_dt_sq;


        r[i].x      += ddx;
        r[i].y      += ddy;
        r[i].z      += ddz;

        v_norm.x = disp1*(v_norm.x + (f_i.x + gamma*(v_hydr.x-v_norm.x))*half_dt);
        v_norm.y = disp1*(v_norm.y + (f_i.y + gamma*(v_hydr.y-v_norm.y))*half_dt);
        v_norm.z = disp1*(v_norm.z + (f_i.z + gamma*(v_hydr.z-v_norm.z))*half_dt);

        v[i].x        = v_norm.x;
        v[i].y        = v_norm.y;
        v[i].z        = v_norm.z;
      }
    return( parameters.mass/2.0);
} 


/* ########################################################################## */
void      verlet_pt2_1filament(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *fa,
                       vec_s  *vhydro,
                       param_s parameters )
       
/* ########################################################################## */
{
  static char  init_flag = TRUE;
  static float dt, half_dt, half_dt_sq, disp1;
  static float gamma;
  static float bond, hydr_radius;
  static vec_s *v_hyd = NULL;
  static int n_mon1, n_mon2, n_fil;
  int  i;

  vec_s f_i, v_norm, v_hydro;

    if( init_flag )
    {
        dt         = parameters.dt;
        bond       = parameters.b_l;
        hydr_radius = parameters.hydr_radius;
        half_dt    = dt/2.0;
        half_dt_sq = dt*dt/2.0;
        n_mon1      = parameters.n_mon1;
	gamma = parameters.gamma1;


        init_flag  = FALSE;
    }
    if(v_hyd == NULL)
      v_hyd = (vec_s *)malloc(n_mon1*sizeof( vec_s ));

    if (parameters.oseen ==1){
      oseen_1filament(r, fa, v_hyd, parameters);
    }
    if (parameters.blake ==1){
      blake_1filament(r, fa, v_hyd, parameters);
    }
  
    for( i = 0; i < n_mon1; i++ )
    {

        disp1         = 1.0/(1.0 + 0.5*gamma*dt);
        f_i.x      = fa[i].x;
        f_i.y      = fa[i].y;
        f_i.z      = fa[i].z;
   

        v_norm.x   = v[i].x;
        v_norm.y   = v[i].y;
        v_norm.z   = v[i].z;
  
        v_hydro.x   = v_hyd[i].x;
        v_hydro.y   = v_hyd[i].y;
        v_hydro.z   = v_hyd[i].z;

        vhydro[i].x   = v_hyd[i].x;
        vhydro[i].y   = v_hyd[i].y;
        vhydro[i].z   = v_hyd[i].z;


        v_norm.x = disp1* (f_i.x + gamma*(v_hydro.x))*half_dt;
        v_norm.y = disp1* (f_i.y + gamma*(v_hydro.y))*half_dt;
        v_norm.z = disp1* (f_i.z + gamma*(v_hydro.z))*half_dt;
       
        v[i].x      += v_norm.x;
        v[i].y      += v_norm.y;
        v[i].z      += v_norm.z;
    }
}



/* ########################################################################## */
void     output_data_1filament(  vec_s  *r,
                       vec_s  *v,
                       vec_s  *f,
           		param_s params,
                       int     step,
                       float   *fper_ave,
                       float   *fpar_ave       )
/* ########################################################################## */
{
    FILE  *out_file_1, *out_file_2, *out_file_3;
    char  file_name[MAX_BUFF];
    int   i,j, n_mon1, n_mon2;
    float rijx, rijy, rijz, rij, t, frac, dprod;
    float f_parx_tot, f_pary_tot, f_parz_tot;
    float f_perx_tot, f_pery_tot, f_perz_tot;
    float hydr_radius, gamma, mag;

    vec_s f_par, f_per;

    
    t    = step*params.dt;
    frac = t*params.freq/(2.0*PI); 
    frac = step/params.n_cycle;
    n_mon1 = params.n_mon1;


    sprintf( file_name,"%s%s%.2f", params.dir_name, "out_", frac );
    out_file_1 = fopen( file_name, "w" );

    for( i = 0; i < n_mon1; i++ )
  {
    fprintf( out_file_1, "%.10e %.10e %.10e %.10e %.10e %.10e \n", 
       r[i].x, r[i].y, r[i].z, v[i].x, v[i].y, v[i].z );
  }

    
    fclose( out_file_1 );

} 




