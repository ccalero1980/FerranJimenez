typedef struct {
int   n_steps;
int   n_mon;
int   n_mon1;
int   n_mon2;
int   n_cycle;
int   n_fil;/*   Mar */
int   blake;
int   oseen;
float mass;
float mass1;
float mass2;
float gamma;
float gamma_head;
float viscosity;
float gamma1;
float gamma2;
float hydr_radius;
float dt;
float k1;
float k;
float phi_0;
float b_l;
float tolerance;
float amp;
float freq;
float w_length;
int   n_left;
int   n_right;
int   out_freq;
char  infile_name[MAX_BUFF];
char  outfile_name[MAX_BUFF];
char  dir_name[MAX_BUFF];
float amplitude; /* Fab */
float Sperm_Number; /* Fab */
float dist; /*  Mar */
float Fx;
float Fz;
float Bx;
float Bz;
float chi;
float radius_para;
float m;
float mu0;
float gravity;
} param_s;


void   configure_sys( int,
                      char    *[],
                      param_s *,
                      vec_s   **,
                      vec_s   **);

void   simulate( param_s,
                 vec_s *,
                 vec_s *  );



void  steric_force( vec_s  , 
                    float  ,
                    vec_s * );

void  steric_force2( vec_s  , 
                    float  ,
                    vec_s * );

/* Bending force routine
Compute bending forces. 
It calls angle_force
Return:   u (???),
          f (forces)
          others?
*/
float  bending_forces_1filament( vec_s *,
                       float *,
                       vec_s *,
                       param_s  );

/**************************
Computes angular forces,
diedral????
The forces are combinations of distances......????
*/
float  angle_force( vec_s,
                    vec_s,
                    float,
                    float,
                    float,
                    float,
                    vec_s * );

void    constrain_positions_1filament(  vec_s *,
                              vec_s *,
                              vec_s *,
                              vec_s *,
                              param_s );

void   constrain_velocities_1filament(  vec_s *, 
                              float *,
                              vec_s *,
                              vec_s *,            
                              param_s );

void   constrain_forces_1filament(      vec_s *, 
                              float *,
                              vec_s *,
                              vec_s *,
                              param_s );

float   verlet_pt1_1filament( vec_s *,
                    vec_s *,
                    vec_s *,
                    vec_s *,
                    param_s );

void   verlet_pt2_1filament( vec_s *,
                   vec_s *,
                   vec_s *,
                   vec_s *,
                   param_s );



void oseen_1filament( vec_s *,
      vec_s *,
      vec_s *,
      param_s);


void blake_1filament( vec_s *,
      vec_s *,
      vec_s *,
      param_s);



void  output_data_1filament( vec_s *,
                   vec_s *,
                   vec_s *,
                   param_s,
                   int,
                   float *,
                   float *);



/* Fab */
void print_force(vec_s *, param_s, int, char *);
