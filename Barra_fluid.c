#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/* #include <ran.h>*/
#include "useful.h"
#include "Funcions/os3many_multiple_proves.h"

#include "Funcions/vectors.h"
#include "Funcions/vectors.c"
#include "Funcions/funcions_proves.c"

#define BFIELD
//#define SQUARE
#define CONSTANT_BFIELD
#define INDUCTION
#define STERIC
//#define GRAVITY
#define HYDROCONSISTENT // Sets the gamma coefficient of the beads composing the filament in a way consistent with the hydrodynamics of a slender body (see Bailey_PRE2009)

//#define BUCKLING_INDUIT
//#define BUCKLING
//#define TRANS
//#define CONSTANT_FORCE
//#define HEAD
int main( int argc, char *argv[] )
{
    vec_s *r, *v;

    param_s parameters;
    configure_sys( argc, argv, &parameters, &r, &v);
    simulate( parameters, r, v );

    return(1);
}

/* ######################################################################## */
void  configure_sys( int      argc,
                     char    *argv[],
                     param_s *parameters,
                     vec_s  **r,
                     vec_s  **v  )
/* ######################################################################### */
{
    FILE  *param_file, *info_file, *config_file;

    char   buffer[MAX_BUFF], file_name[MAX_BUFF], dir_name[MAX_BUFF];
    float  sp_numb4, tinercia, Fmax;
    int    n_mon1, n_fil,i, i0;
    float  rijx, rijy, rijz, B;

    seed = -789219;

    if( argc == 1 )
    {
              //printf("\nNo parameter file name\n" );
       exit( 1 );
    }
    param_file = fopen( argv[1], "r" );
    if( param_file == NULL )
    {
         //printf("\nParameter file not found\n");
         exit( 1 );
    }
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%s", file_name );
    config_file = fopen( file_name, "r" );
    if( config_file == NULL )
    {
         //printf("\nConfiguration file %s not found\n", file_name );
         fflush( stdout );
         exit( 1 );
    }

    sprintf( parameters->outfile_name,"%s%s", file_name, "_new" );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%s", dir_name );
    sprintf( file_name,"%s%s", dir_name, "info");
    sprintf( parameters->dir_name,"%s", dir_name );
    info_file = fopen( file_name, "w" );
    if( info_file == NULL )
    {
         //printf("\nCannot access directory for output %s\n", file_name );
         fflush( stdout );
         exit( 1 );
    }
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->n_steps) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->dt) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->k1) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->mass1) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->viscosity) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->gamma1) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->hydr_radius) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->tolerance) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->n_cycle) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->Bx) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->Bz) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->freq) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->m) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->mu0) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%le", &(parameters->gravity) );
    fgets( buffer, MAX_BUFF, param_file );

    sscanf( buffer, "%d", &(parameters->oseen) );
    fgets( buffer, MAX_BUFF, param_file );
    sscanf( buffer, "%d", &(parameters->blake) );

    printf("parameters in\n");
   

    n_mon1 = 0;

    while( fgets( buffer, MAX_BUFF, config_file ) != NULL )
   {
         n_mon1++;
    }
    n_mon1 = n_mon1;
    
    rewind(config_file);

    fflush( stdout );
    
    /*   posiz, velocita' di tutti i bead */  

    (*r) = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
    (*v) = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
      for(i=0;i<n_mon1;i++){     
        fgets(  buffer, MAX_BUFF, config_file );
        sscanf( buffer,"%le%le%le%le%le%le", &((*r)[i].x),
                                                &((*r)[i].y),
                                                &((*r)[i].z),
                                                &((*v)[i].x),
                                                &((*v)[i].y),
                                                &((*v)[i].z));    
                  //printf("%le\t%le\n", (*r)[i].x, (*v)[i].x);

    
      }
    fclose( config_file );

    
    // //printf("n_mon %d, n_fil %d\n",n_mon, n_fil);

    rijx = (*r)[0].x - (*r)[1].x;
    rijy = (*r)[0].y - (*r)[1].y;
    rijz = (*r)[0].z - (*r)[1].z;
    parameters->b_l         = sqrt(rijx*rijx + rijy*rijy + rijz*rijz );
    parameters->n_mon1       = n_mon1;
    parameters->mass1       /= (float)n_mon1;



    parameters->gamma1      /= parameters->mass1*n_mon1;
#ifdef HYDROCONSISTENT
    parameters->gamma1      = 6.*PI*parameters->viscosity*parameters->b_l*0.559/parameters->mass1;
#endif
    parameters->k1          /= parameters->b_l*parameters->b_l;


        
    fprintf( info_file, "Number of steps   = %d\n", parameters->n_steps  );
    fprintf( info_file, "Time-step         = %e\n", parameters->dt       );
    fprintf( info_file, "Number of monomers   = %d\n", parameters->n_mon1  );
    fprintf( info_file, "Bending modulus   = %e\n", parameters->k1        );
    fprintf( info_file, "Bond length       = %e\n", parameters->b_l      );
    fprintf( info_file, "Tolerance         = %e\n", parameters->tolerance); 
    fprintf( info_file, "Friction gamma    = %e\n", parameters->gamma1    );
    fprintf( info_file, "Frequency         = %e\n", parameters->freq*parameters->dt/(2.*PI)    );
    fprintf( info_file, "Period         = %e\n\n (steps)", 2*PI/(parameters->freq*parameters->dt));
 
    fprintf( info_file, "mu 0 = %e\n", parameters->mu0); 



    if (parameters->oseen == 1){
      fprintf( info_file, "Oseen \n");
    }
    if (parameters->blake == 1){
      fprintf( info_file, "Blake \n");
    }
    fclose(info_file);
}

/* ########################################################################## */

void  simulate(     param_s  params,
                     vec_s  *r,
                     vec_s  *v)
/* ########################################################################## */
{
    vec_s *fc_1, *fc_2, *fc_3, *v_sv;
    vec_s *f_bend, *fc_1a, *fc_2a,*fc_3a, *f_dip, *f_steric, *f_torque;
    vec_s *vhydro;


    int    i, j, k, l,m, n_mon1, step, n_fil, nmig;

    float  *mags, *mags1;
    float  fx_tot, fy_tot, fz_tot, u, ke, kk;
    float  rijx, rijy, rijz;
    float  fpar_ave, fper_ave;
    float  cycle_time, H, shear_rate, Torque_Poiseuille;
    vec_s  r0f, m1, m2, Bfield, Bin, Bin0;
    vec_s *rij, *rij1, *f,*fa, *tang,*v1, *r1, *rm, *rcm;

    FILE *conf = fopen("conf", "w");
   // FILE *flog_fa;
    /* Fab */ float freq;


    /* MCL */ float  b_l ;
    /* MCL */ float  dr0x, dr0y, dr0z,drijx,drijy, drijz;
    /* MCL */ float  Torque, amin, zmin;

    FILE *frcm = fopen("RCM.dat", "w");
    FILE *flog = fopen("fil.log", "w");

 
      
    /* Fab */ freq = params.freq;


    
    n_mon1 = params.n_mon1;   
    
/* Forze  a = di tutti  */    

    f_bend   = (vec_s *)malloc( n_mon1*sizeof( vec_s ));
    f_dip   = (vec_s *)malloc(sizeof( vec_s ));
    f_steric   = (vec_s *)malloc(sizeof( vec_s ));
    f_torque   = (vec_s *)malloc(sizeof( vec_s ));

    fc_1a    = (vec_s *)malloc( n_mon1*sizeof( vec_s ));  
    fc_1     = (vec_s *)malloc( n_mon1*sizeof( vec_s ));

    fc_2a    = (vec_s *)malloc( n_mon1*sizeof( vec_s ));
    fc_2     = (vec_s *)malloc( n_mon1*sizeof( vec_s ));
 
    fc_3a    = (vec_s *)malloc( n_mon1*sizeof( vec_s ));
    fc_3     = (vec_s *)malloc( n_mon1*sizeof( vec_s ));

    v_sv     = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
    vhydro     = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
    
    

    rij   = (vec_s *)malloc((n_mon1+2)*sizeof( vec_s ));
    rij1  = (vec_s *)malloc((n_mon1 + 1)*sizeof( vec_s ));

    fa    = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
    f     = (vec_s *)malloc(n_mon1*sizeof( vec_s ));

    v1    = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
    r1    = (vec_s *)malloc(n_mon1*sizeof( vec_s ));

    mags  = (float *)malloc(n_mon1*sizeof( float ));
    mags1 = (float *)malloc(n_mon1*sizeof( float ));
    tang  = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
  
    rm    = (vec_s *)malloc(n_mon1*sizeof( vec_s ));
    rcm = (vec_s *)malloc(sizeof( vec_s ));

   printf(" inici programa \n");

    
    cycle_time  = params.n_cycle*params.dt;

    for( i = 1; i < n_mon1; i++ )
    {
       rijx     = r[i-1].x - r[i].x;
       rijy     = r[i-1].y - r[i].y;
       rijz     = r[i-1].z - r[i].z;
      if(i%n_mon1==0){
      	 rijx     = 0;
      	 rijy     = 0;
      	 rijz     = 0;  
    	} 
     
      mags[i]  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
      rij[i].x = rijx;
      rij[i].y = rijy;
      rij[i].z = rijz;


    }


    for( i = 0; i < n_mon1; i++ )
    {
       fa[i].x = 0.0;
       fa[i].y = 0.0;
       fa[i].z = 0.0;

       fc_1a[i].x = 0.0;
       fc_1a[i].y = 0.0;
       fc_1a[i].z = 0.0;
       fc_2a[i].x = 0.0;
       fc_2a[i].y = 0.0;
       fc_2a[i].z = 0.0;
       fc_3a[i].x = 0.0;
       fc_3a[i].y = 0.0;
       fc_3a[i].z = 0.0;
    }


    for( i = 0; i < n_mon1; i++ ){
    	rij1[i].x = rij[i].x;
    	rij1[i].y = rij[i].y;
    	rij1[i].z = rij[i].z;


    	mags1[i] = mags[i];
    	
    	f[i].x = fa[i].x;
    	f[i].y = fa[i].y;
    	f[i].z = fa[i].z;
    }

    u    = bending_forces_1filament( rij1, mags1, f, params );
     
    for( i = 0; i < n_mon1; i++ ){
    	fa[i].x += f[i].x;
    	fa[i].y += f[i].y;
    	fa[i].z += f[i].z;

	
    }

#ifdef GRAVITY
    for( i = 0; i < n_mon1; i++ ){
        fa[i].z -= params.gravity; // gravity g=10.m/s^2 in units of Lfilament/t0^2 
    }
#endif

#ifdef STERIC
    zmin = params.b_l; // minimum distance of filament with wall at z = 0

//    STERIC INTERACTION WITH WALL AT Z=0 -----------------
    for( i = 0; i < n_mon1; i++ ){
        f_steric[0].x = 0.0;
        f_steric[0].y = 0.0;
        f_steric[0].z = 0.0;

        r0f.x = 0.0;
        r0f.y = 0.0;
        r0f.z = r[i].z;


        steric_force2(r0f, zmin, f_steric);
        fa[i].x += f_steric[0].x/params.mass1;
        fa[i].y += f_steric[0].y/params.mass1;
        fa[i].z += f_steric[0].z/params.mass1;

    }

#endif

#ifdef BFIELD

      // BFIELD
      Bfield.x = params.Bx*sin(freq*0.*params.dt);
      Bfield.y = 0.0;
      Bfield.z = params.Bz*cos(freq*0*params.dt);
      //   Bfield.z = 0.0;
#ifdef SQUARE
      Bfield.x = params.Bx;
      Bfield.y = 0.0;
      Bfield.z = 0.0;
#endif

#ifdef CONSTANT_BFIELD
      Bfield.x = params.Bx;
      Bfield.y = 0.0;
      Bfield.z = params.Bz;
#endif
Bin.x = 0.0;
Bin.y = 0.0;
Bin.z = 0.0;


    // INTERACTION WITH BFIELD   
      for( i = 1; i < n_mon1; i++ ){
          f_torque[0].x = 0.0;
          f_torque[0].y = 0.0;
          f_torque[0].z = 0.0;

          fa[i].x += f_torque[0].x;
          fa[i].y += f_torque[0].y;
          fa[i].z += f_torque[0].z;

        } 


#endif


////////////////////////////////////////////////////////////////////////////////////////////
//COMENÇA EL CICLE TEMPORAL!!!!!!!!!!!!!!!
////////////////////////////////////////////////////////////////////////////////////////////

  for( j = 0; j < params.n_steps; j++ ){
    printf("step j = %i \n",j);
   

    for( i = 0; i < n_mon1; i++ ){
  	  rm[i].x=0.0;
  	  rm[i].y=0.0;
  	  rm[i].z=0.0;
    }


    rcm[0].x=0.0;
    rcm[0].y=0.0;
    rcm[0].z=0.0;


    for( k = 0; k < params.n_cycle; k++ ){
      step = j*params.n_cycle + k;
	  

      //fflush( stdout );


	    for( i = 0; i < n_mon1; i++ ){
	      rij1[i].x = rij[i].x;
	      rij1[i].y = rij[i].y;
	      rij1[i].z = rij[i].z;
	 
	   
	      mags1[i] = mags[i];
	    
	      v1[i].x = v[i].x;
	      v1[i].y = v[i].y;
	      v1[i].z = v[i].z;
	    }
	     constrain_velocities_1filament( rij1, mags1, v1,fc_1, params );
	
	    for( i = 0; i < n_mon1; i++ ){
	      v[i].x = v1[i].x;
	      v[i].y = v1[i].y;
	      v[i].z = v1[i].z;


	      fc_1a[i].x   = fc_1[i].x;
	      fc_1a[i].y   = fc_1[i].y;
	      fc_1a[i].z   = fc_1[i].z;
	    }
	  
      ke = verlet_pt1_1filament( r, v, fa, vhydro, params );

	    for( i = 0; i < n_mon1; i++ ){

	      rij1[i].x = rij[i].x;
	      rij1[i].y = rij[i].y;
	      rij1[i].z = rij[i].z;
	      
	      mags1[i] = mags[i];

	      v1[i].x = v[i].x;
	      v1[i].y = v[i].y;
	      v1[i].z = v[i].z;

	      r1[i].x = r[i].x;
	      r1[i].y = r[i].y;
	      r1[i].z = r[i].z;
	    }    
	    constrain_positions_1filament( r1, rij1, v1, fc_2, params );

	    
	    for( i = 0; i < n_mon1; i++ ){

  	    v[i].x = v1[i].x;
  	    v[i].y = v1[i].y;
  	    v[i].z = v1[i].z;
  	    
  	   
  	    r[i].x = r1[i].x;
  	    r[i].y = r1[i].y;
  	    r[i].z = r1[i].z;

  	    fc_2a[i].x   = fc_2[i].x;
  	    fc_2a[i].y   = fc_2[i].y;
  	    fc_2a[i].z   = fc_2[i].z;
	    }

	    for( i = 1; i < n_mon1; i++ ){
	      rijx     = r[i-1].x - r[i].x;
	      rijy     = r[i-1].y - r[i].y;
	      rijz     = r[i-1].z - r[i].z;
	      if(i%n_mon1==0){
      		rijx     = 0;
      		rijy     = 0;
      		rijz     = 0;
	      }
	      mags[i]  = sqrt(rijx*rijx + rijy*rijy + rijz*rijz);
	      rij[i].x = rijx;
	      rij[i].y = rijy;
	      rij[i].z = rijz;
	    }


	
      for( i = 0; i < n_mon1; i++ ){
    	    fa[i].x = 0.0;
    	    fa[i].y = 0.0;
    	    fa[i].z = 0.0;
    	}
#ifdef GRAVITY
    for( i = 0; i < n_mon1; i++ ){
        fa[i].z -= params.gravity; // gravity g=10.m/s^2 in units of Lfilament/t0^2 
    }
#endif

#ifdef STERIC
      zmin = params.b_l; // minimum distance of filament with wall at z = 0

/*
//    STERIC INTERACTION WITH WALL AT Z=0 -----------------
      for( i = 0; i < n_mon1; i++ ){
          f_steric[0].x = 0.0;
          f_steric[0].y = 0.0;
          f_steric[0].z = 0.0;

          r0f.x = 0.0;
          r0f.y = 0.0;
          r0f.z = r[i].z;

          steric_force2(r0f, zmin, f_steric);
          fa[i].x += f_steric[0].x/params.mass1;
          fa[i].y += f_steric[0].y/params.mass1;
          fa[i].z += f_steric[0].z/params.mass1;

      }
*/

// ------------------------------------------------------------------

#endif 

#ifdef BFIELD

      // BFIELD
      Bfield.x = params.Bx*sin(freq*step*params.dt);
    //  Bfield.x = params.Bx*cos(freq*step*params.dt);
      Bfield.y = 0.0;
      Bfield.z = params.Bz*cos(freq*step*params.dt);
    //  Bfield.z = params.Bz*sin(freq*step*params.dt);
      //   Bfield.z = 0.0;
#ifdef SQUARE
      if(sin(freq*step*params.dt) > 0.0){
        Bfield.x = params.Bx;
      }
      else if(sin(freq*step*params.dt) < 0.0){
        Bfield.x = -params.Bx;
      }
      if(cos(freq*step*params.dt) > 0.0){
        Bfield.z = params.Bz;
      }
      else if(cos(freq*step*params.dt) < 0.0){
        Bfield.z = -params.Bz;
      }
      Bfield.y = 0.0;
#endif

#ifdef CONSTANT_BFIELD
      Bfield.x = params.Bx;
      Bfield.y = 0.0;
      Bfield.z = params.Bz;
#endif

 

    // INTERACTION BY BFIELD ON THE MAGNETIC BAR  

    // INCLOURE !!!


      for( i = 1; i < n_mon1; i++ ){
          f_torque[0].x = 0.0;
          f_torque[0].y = 0.0;
          f_torque[0].z = 0.0;

          fa[i].x += f_torque[0].x;
          fa[i].y += f_torque[0].y;
          fa[i].z += f_torque[0].z;

        } 


#endif



	    for( i = 0; i < n_mon1; i++ ){      
	      f_bend[i].x = 0.0;
	      f_bend[i].y = 0.0;
	      f_bend[i].z = 0.0;
	
	      rij1[i].x = rij[i].x;
	      rij1[i].y = rij[i].y;
	      rij1[i].z = rij[i].z;

	      mags1[i] = mags[i];
	    }

	    u  = bending_forces_1filament( rij1, mags1, f_bend, params );
     
	    for( i = 0; i < n_mon1; i++ ){
	      fa[i].x += f_bend[i].x;
	      fa[i].y += f_bend[i].y;
	      fa[i].z += f_bend[i].z;	      
	    }
	  

	    for( i = 0; i < n_mon1; i++ ){
	      rij1[i].x = rij[i].x;
	      rij1[i].y = rij[i].y;
	      rij1[i].z = rij[i].z;
	 	   
	      mags1[i] = mags[i];
	    
	      f[i].x = fa[i].x;
	      f[i].y = fa[i].y;
	      f[i].z = fa[i].z;	
	    }
        

	    constrain_forces_1filament(rij1, mags1, f, fc_3, params );
	
	    for( i = 0; i < n_mon1; i++ ){
	      fc_3a[i].x   = fc_3[i].x;
	      fc_3a[i].y   = fc_3[i].y;
	      fc_3a[i].z   = fc_3[i].z;
	    }
	  
        for( i = 0; i < n_mon1; i++ )
          {

	    fa[i].x  += fc_3a[i].x;
	    fa[i].y  += fc_3a[i].y;
	    fa[i].z  += fc_3a[i].z;  

        }       
        verlet_pt2_1filament( r, v, fa , vhydro, params );




	  /*Position of the CM*/
	

	for( i = 0; i < n_mon1; i++ ){
	  rm[i].x += r[i].x;
	  rm[i].y += r[i].y;
	  rm[i].z += r[i].z;
    }


    if(step%params.n_cycle==0){
          rcm[0].x = 0.0;
          rcm[0].y = 0.0;
          rcm[0].z = 0.0;

   
        for(m=0;m<n_mon1; m++){
          rcm[0].x += r[m].x/n_mon1;
          rcm[0].y+= r[m].y/n_mon1;
          rcm[0].z += r[m].z/n_mon1;
        }


	  fprintf(frcm,"%i\t%le\t%le\t%le\t%le\t%le\t%le\n", step, rcm[0].x, rcm[0].y,rcm[0].z, Bfield.x, Bfield.y, Bfield.z );


    }


       
       output_data_1filament( r, v, fa, params, step, &fper_ave, &fpar_ave );
       fflush( stdout );
	  



   }
  }
fclose(frcm);
fclose(flog);


}


