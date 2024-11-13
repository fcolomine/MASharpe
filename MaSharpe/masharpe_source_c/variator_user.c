/*========================================================================
 MASHARPE
     
  ========================================================================
*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include "variator.h"
#include "variator_user.h"

/*--------------------| global variable definitions |-------------------*/

/*==== declared in variator_user.h used in other files as well ====*/

char *log_file = "dtlz_error.log"; /**** Changed for DTLZ.*/

char paramfile[FILE_NAME_LENGTH]; /* file with local parameters */

char *nro_corridas = "corridas.log";


//char modelo[FILE_NAME_LENGTH]; 

/*==== only used in this file ====*/

/* local parameters from paramfile*/
char problem[FILE_NAME_LENGTH]; 
int seed;   /* seed for random number generator */
int number_decision_variables; /* length of the binary string */
int maxgen; /* maximum number of generations (stop criterion) */
int gen=0;
int impgen=0;
char outfile[FILE_NAME_LENGTH]; /* output file for last population */
double individual_mutation_probability;
double individual_recombination_probability;
double variable_mutation_probability;
double variable_swap_probability;
double variable_recombination_probability;
double eta_mutation;
double eta_recombination;
int use_symmetric_recombination;
int Card;
double randcard;
int corridas;
double nrocorridas=0.0;
int nrocorridasBL=0;
int nrocorridasBLV=0;
int nrocorridasmax=0;
double hill;
int iteracion;
double hill1;
double hill2;
double hill3;
int iteracion1;
int iteracion2;
int diversidad;
int diversidaddinamica;
ProblemInfo* p;
double tprom,maximosharpe,hipervolumemax,hipervolumemin;
//double *minimo,*maximo,*prom;
double* minimo;
double* maximo;
double* prom;
double* maxxsharpe;
double mejores[50][30];

/*-------------------------| individual |-------------------------------*/

void free_individual(individual *ind) 
/* Frees the memory for one indiviual.

   post: memory for ind is freed
*/
{
      
     /**********| added for DTLZ |**************/
    
     if (ind == NULL)
          return;
     
     free(ind->x);
     free(ind->f);
     free(ind->h);
     free(ind->Phenotype);
     /**********| addition for DTLZ end |*******/
     
     free(ind);
}

double get_objective_value(int identity, int i)
/* Gets objective value of an individual.

   pre: 0 <= i <= dimension - 1 (dimension is the number of objectives)

   post: Return value == the objective value number 'i' in individual '*ind'.
         If no individual with ID 'identity' return value == -1. 
*/   
{
    
     /**********| added for DTLZ |**************/
     individual *temp_ind;
     /**********| addition for DTLZ end |*******/
     
     double objective_value = -1.0;

     assert(0 <= i && i < dimension); /* asserting the pre-condition */
     
     /**********| added for DTLZ |**************/
    
     if (i < 0 || i > (dimension - 1))
	  return(-1);
     
     temp_ind = get_individual(identity);
     if (temp_ind == NULL)
	 return(-1);
     
     objective_value = temp_ind->f[i];
     
     /**********| addition for DTLZ end |*******/
  
     return (objective_value);
}

/*-------------------------| statemachine functions |-------------------*/

int state0() 
/* Do what needs to be done in state 0.

   pre: The global variable 'paramfile' contains the name of the
        parameter file specified on the commandline.
        The global variable 'alpha' contains the number of indiviuals
        you need to generate for the initial population.
                
   post: Optionally read parameter specific for the module.
         Optionally do some initialization.
         Initial population created.
         Information about initial population written to the ini file
         using write_ini().
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     /**********| added for DTLZ |**************/
     int i;
     /**********| addition for DTLZ end |*******/
    int *parent_identities1;
    parent_identities1 = (int *) malloc(alpha * sizeof(int));
    p = (ProblemInfo*)malloc(sizeof(ProblemInfo));
     int result; /* stores return values of called functions */
     int *initial_population; /* storing the IDs of the individuals */
     initial_population = (int *) malloc(alpha * sizeof(int)); 
     if (initial_population == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }

     /**********| added for DTLZ |**************/
     result = read_local_parameters();
/*--------------------------AGREGADO POR GASSHARPE---------*/  
   
     Leerfondos(p,problem);
/*--------------------------AGREGADO POR GASSHARPE---------*/ 
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't read local parameters");
          return (1);
     }
     
     /* initializing the first alpha individuals */
     for(i = 0; i < alpha; i++)
     {
     
 	 individual *ind = new_individual();
 	 eval(ind);
	 initial_population[i] = add_individual(ind);
         if(initial_population[i] == -1)
            return(1);
     } 
     if ((drand(1)<hill2)&&(nrocorridas<(double)nrocorridasmax)){
        printf("%s","PASO POR EL ESTADO CORRECTOR");
        gen=0;
        result = variate1(initial_population,parent_identities1);
        printf("%s","SALIO DEL ESTADO CORRECTOR");
        //getch();
        if (result != 0) return (1);
        
        result = write_ini(parent_identities1);
          if (result != 0)
                     {                
                      log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write ini");
                      free(initial_population);
                       return (1);
                       }
          free(parent_identities1);   
     }else{
           result = write_ini(initial_population);
           if (result != 0)
              { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write ini");
          free(initial_population);
          return (1);
          }
          
          free(initial_population);
          }; 
     gen = 1;
     printf("%s","PASO POR EL ESTADO 0");
     /**********| addition for DTLZ end |*******/

     /*result = write_ini(initial_population);
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write ini");
          free(initial_population);
          return (1);
     }

     free(initial_population);*/
     //getch();
     wait(0.1);
     return (0);
}



int state2()
/* Do what needs to be done in state 2.

   pre: The global variable 'mu' contains the number of indiviuals
        you need to read using 'read_sel()'.
        The global variable 'lambda' contains the number of individuals
        you need to create by variation of the individuals specified the
        'sel' file.
        
   post: Optionally call read_arc() in order to delete old uncessary
         individuals from the global population.
         read_sel() called
         'lambda' children generated from the 'mu' parents
         Children added to the global population using add_individual().
         Information about children written to the 'var' file using
         write_var().
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
      
     int *parent_identities, *offspring_identities; /* array of identities */
     int result; /* stores return values of called functions */

     parent_identities = (int *) malloc(mu * sizeof(int)); 
     if (parent_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }

     offspring_identities = (int *) malloc(lambda * sizeof(int)); 
     if (offspring_identities == NULL)
     {
          log_to_file(log_file, __FILE__, __LINE__, "variator out of memory");
          return (1);
     }
     
     result = read_sel(parent_identities);
     if (result != 0) /* if some file reading error occurs, return 2 */
          return (2);

     result = read_arc(); 
    if (result != 0) /* if some file reading error occurs, return 2 */
          return (2);

     /**********| added for DTLZ |**************/
     
     

     if (((drand(1)<hill||hill==0.000123)&&((gen>iteracion)&&(gen<iteracion2))&&(nrocorridas<(double)nrocorridasmax))){
        if (!is_finished()){                                                                                    
        result = variate1(parent_identities, offspring_identities);
        printf("%s","VARIADOR HIBRIDO...."); 
        if (result != 0)
            return (1);
        };}
     else{
          if (!is_finished()){     
          result = variate(parent_identities, offspring_identities);
          printf("%s","VARIADOR NATURAL...."); 
          if (result != 0)
             return (1);   
            }; };
     if(hill>0){
                if(nrocorridas<(double)nrocorridasmax){
                 gen++;
                 }else{
                       printf("Me pare en la generacion %i",gen);
                       //getch(); 
                       impgen=gen;
                       gen=maxgen;};
     }else{
                       gen++;
                       };
                      
     /**********| addition for DTLZ end |*******/


     result = write_var(offspring_identities);
     
     if (result != 0)
     { 
          log_to_file(log_file, __FILE__, __LINE__,
                      "couldn't write var");
          free(offspring_identities);
          free(parent_identities);
          return (1);
     }

     free(offspring_identities);
     free(parent_identities);
     
     return (0);
}
 

int state4() 
/* Do what needs to be done in state 4.

   pre: State 4 means the variator has to terminate.

   post: Free all memory.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     /**********| added for DTLZ |**************/
     
     int result;
     result = read_arc();

     if (0 == result) /* arc file correctly read
                         this means it was not read before,
                         e.g., in a reset. */
     {
        write_output_file();
     }
     
     /**********| addition for DTLZ end |*******/
     
     return (0);
}


int state7()
/* Do what needs to be done in state 7.

   pre: State 7 means that the selector has just terminated.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     return(0);  
}


int state8()
/* Do what needs to be done in state 8.

   pre: State 8 means that the variator needs to reset and get ready to
        start again in state 0.

   post: Get ready to start again in state 0. 
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     /**********| added for DTLZ |**************/

   int result;
   
   gen = 1;
     
   result = read_arc();

   if (0 == result) /* arc file correctly read
                       this means it was not read before */
   {
      write_output_file();
   }
   
     /**********| addition for DTLZ end |*******/
     
     return (0);
}


int state11()
/* Do what needs to be done in state 11.

   pre: State 11 means that the selector has just reset and is ready
        to start again in state 1.

   post: You probably don't need to do anything, just return 0.
         Return value == 0 if successful,
                      == 1 if unspecified errors happened,
                      == 2 if file reading failed.
*/
{
     return (0);  
}


int is_finished()
/* Tests if ending criterion of your algorithm applies.

   post: return value == 1 if optimization should stop
         return value == 0 if optimization should continue
*/
{
         
     /**********| added for DTLZ |**************/
     return (gen >= maxgen);
     /**********| addition for DTLZ end |*******/
}


/**********| added for DTLZ |**************/

int read_local_parameters()
{
     FILE *fp;
     char str[CFG_NAME_LENGTH];

     /* reading parameter file with parameters for selection */
     fp = fopen(paramfile, "r"); 
     assert(fp != NULL);

     if(dimension < 0)
     {
          log_to_file(log_file, __FILE__, 
                      __LINE__, "can't handle that dimension");
          return(1);
     } 

     if(mu != lambda)
     {
          log_to_file(log_file, __FILE__, 
                      __LINE__, "can't handle mu != lambda");
          return(1);
     }

     fscanf(fp, "%s", str);
     assert(strcmp(str, "problem") == 0);
     fscanf(fp, "%s", problem); /* fscanf() returns EOF if
                                   reading failed. */
      //sprintf(problem, "%s", modelo);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "seed") == 0);
     fscanf(fp, "%d", &seed);
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "number_decision_variables") == 0);
     fscanf(fp, "%d", &number_decision_variables);
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "maxgen") == 0);
     fscanf(fp, "%d", &maxgen);
     maxgen=gene;
     fscanf(fp, "%s", str);
     assert(strcmp(str, "outputfile") == 0);
     fscanf(fp, "%s", outfile); /* fscanf() returns EOF if
                                   reading failed. */

     fscanf(fp, "%s", str);
     assert(strcmp(str, "individual_mutation_probability") == 0);
     fscanf(fp, "%le", &individual_mutation_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "individual_recombination_probability") == 0);
     fscanf(fp, "%le", &individual_recombination_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "variable_mutation_probability") == 0);
     fscanf(fp, "%le", &variable_mutation_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "variable_swap_probability") == 0);
     fscanf(fp, "%le", &variable_swap_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "variable_recombination_probability") == 0);
     fscanf(fp, "%le", &variable_recombination_probability);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "eta_mutation") == 0);
     fscanf(fp, "%le", &eta_mutation);
     assert(eta_mutation >= 0);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "eta_recombination") == 0);
     fscanf(fp, "%le", &eta_recombination);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "use_symmetric_recombination") == 0);
     fscanf(fp, "%d", &use_symmetric_recombination);

     fscanf(fp, "%s", str);
     assert(strcmp(str, "K") == 0);
     fscanf(fp, "%d", &Card);
     Card=k2;
     fscanf(fp, "%s", str);
     assert(strcmp(str, "random") == 0);
     fscanf(fp, "%le", &randcard);
     assert(randcard >= 0);
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "corridas") == 0);
     fscanf(fp, "%d", &corridas);
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "P_ls") == 0);
     fscanf(fp, "%le", &hill);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "numero_generacion_inicio") == 0);
     fscanf(fp, "%d", &iteracion);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "numero_generacion_final") == 0);
     fscanf(fp, "%d", &iteracion2);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "P_em") == 0);
     fscanf(fp, "%le", &hill1);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "numero_max_itera") == 0);
     fscanf(fp, "%d", &iteracion1);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "numero_de_corridas_max") == 0);
     fscanf(fp, "%d", &nrocorridasmax);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "corrector_pob_inicial") == 0);
     fscanf(fp, "%le", &hill2);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "individual_search") == 0);
     fscanf(fp, "%le", &hill3);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "diversidad") == 0);
     fscanf(fp, "%d", &diversidad);
     printf("%s","Termino de leer");
     srand(seed); /* seeding random number generator */
     
     fclose(fp);

     return (0);
}


/* Performs variation. */
int variate(int *selected, int *result_ids)
{
     int result, i, k;

     result = 1;

     /* copying all individuals from selected */
     for(i = 0; i < mu; i++)
     {
          result_ids[i] = 
               add_individual(copy_individual(get_individual(selected[i])));
          
          if(result_ids[i] == -1)
          {
               log_to_file(log_file, __FILE__, __LINE__,
                           "copying + adding failed");
               return (1);
          }
     }
 
     /* if odd number of individuals, last one is
        left as is */
     if((((double)mu/2) - (int)(mu/2)) != 0) k = mu - 1; 
     else k = mu;
// IMPRIME GENERACION.............. 
int j;
double g=0; 
FILE *fp_out;
char outfile1[FILE_NAME_LENGTH];
individual *temp;
temp=get_individual(result_ids[0]);
sprintf(outfile1, "generacion.log");
fp_out = fopen(outfile1, "a");
assert(fp_out != NULL);
for(i=0;i<mu;i++){
                temp= get_individual(result_ids[i]);  
                fprintf(fp_out, "%d ", i); 
	  for (j = 0; j < dimension; j++)
	  {
	      if (j==0) {
            if(objetivo==1) g=1-(temp->f[j]);
            else g=1000-(temp->f[j]);
            fprintf(fp_out, "%f ", g);}
          else {
               fprintf(fp_out, "%f ", temp->f[j]);
               };
	  };
     fprintf(fp_out, "%1.5f ", (g-p->rf)/sqrt(temp->f[1]));     
	 //printf("%i",gen);
     fprintf(fp_out, "%i ", gen);
     fprintf(fp_out, "%i ", 1);
     fprintf(fp_out, "\n");
   	  
   };
 fclose(fp_out);
 

     /* do recombination */
     for(i = 0; i < k; i+= 2)
     {  
          if (drand(1) <= individual_recombination_probability)
          {
               if (variable_swap_probability > 0)
               {
                    result = uniform_crossover
                      (get_individual(result_ids[i]),
                       get_individual(result_ids[i + 1]));
                    if (result != 0)
                         log_to_file(log_file, __FILE__, 
                                     __LINE__, "recombination failed!");
               }

               if (variable_recombination_probability > 0)
               {
                    result = sbx
                      (get_individual(result_ids[i]), 
                       get_individual(result_ids[i + 1]));
                    if (result != 0)
                         log_to_file(log_file, __FILE__, 
                                     __LINE__, "recombination failed!");
               }
          }
     }
     
     /* do mutation */
     for(i = 0; i < mu; i++)
     {
          if (drand(1) <= individual_mutation_probability)
          { 
               if (variable_mutation_probability > 0)
               {
                    result = mutation(get_individual(result_ids[i]));
                    if(result != 0)
                      log_to_file(log_file, __FILE__, __LINE__,
                                  "mutation failed!");
               }
          }
     }
     
     
           
     
     


  /* do evaluation */
     for(i = 0; i < mu; i++)
     {
          int result;
          result = eval(get_individual(result_ids[i]));
          
     }
    
     
     return (0);
}



/* Performs variation1.*/ 
int variate1(int *selected, int *result_ids)
{
     int result, i, k, n, k1, r;
     result = 1;
double burbuja[50];
    
     for(i = 0; i < mu; i++)
     {
          result_ids[i] = 
               add_individual(copy_individual(get_individual(selected[i])));
          if(result_ids[i] == -1)
          {
               log_to_file(log_file, __FILE__, __LINE__,
                           "copying + adding failed");
               return (1);
          }
     }
 
    
     if((((double)mu/2) - (int)(mu/2)) != 0) k = mu - 1; 
     else k = mu;

int j;
double g=0; 
FILE *fp_out;
char outfile1[FILE_NAME_LENGTH];
individual *temp;
temp=get_individual(result_ids[0]);
sprintf(outfile1, "generacion.log");
fp_out = fopen(outfile1, "a");
assert(fp_out != NULL);
for(i=0;i<mu;i++){
                temp= get_individual(result_ids[i]);  
                fprintf(fp_out, "%d ", i); 
	  for (j = 0; j < dimension; j++)
	  {
	      if (j==0) {
            if(objetivo==1) g=1-(temp->f[j]);
            else g=1000-(temp->f[j]);
            fprintf(fp_out, "%f ", g);}
          else {
               fprintf(fp_out, "%f ", temp->f[j]);
               };
	  };
     fprintf(fp_out, "%1.5f ", (g-p->rf)/sqrt(temp->f[1]));     
	 //printf("%i",gen);
     fprintf(fp_out, "%i ", gen);
     fprintf(fp_out, "%i ", 2);
     fprintf(fp_out, "\n");
   	  
   };
 fclose(fp_out);
 
 //Decision Making 


individual *tmp; 

double rsharpe,sumsharpe,suma,g1,hipervolumemax,hipervolumemin;
minimo = (double*) malloc (sizeof(double)*(p->numfondos));
maximo = (double*) malloc (sizeof(double)*(p->numfondos));
prom = (double*) malloc (sizeof(double)*(p->numfondos));
maxxsharpe = (double*) malloc (sizeof(double)*(p->numfondos));
for(i=0;i<p->numfondos;i++){
                  prom[i]=0;
                  minimo[i]=999;
                  maximo[i]=0;
                  maxxsharpe[i]=0;
              };

for(i=0;i<50;i++){
         for(j=0;j<30;j++){
               
                  mejores[i][j]=0;
                  };
              };
 
int numsharpe,nfondo;
numsharpe=0;
rsharpe=0;
maximosharpe=0;
hipervolumemax=0;
hipervolumemin=999;

tmp=get_individual(result_ids[0]);




            
              
for(i = 0; i < mu; i++)
     {
          suma=0;
          tmp=get_individual(result_ids[i]);
          rsharpe=0;
          //printf("Rf: %1.5f",p->rf);
          if(objetivo==1) g1=1-(tmp->f[0]);
          else g1=1000-(tmp->f[0]);
          //printf("funcion Objetivo 1: %1.5f",g1);
          if(tmp->f[1]>0) rsharpe=(g1-p->rf)/sqrt(tmp->f[1]);
          //printf("rsharpe: %1.5f",rsharpe);
          //getch();
          if(rsharpe>0){
               r=insertarmejores(tmp,rsharpe);
           
                
        
              
          }; 
     };
 //printf("  LLEGUE NUMSHARPE %d",numsharpe); 
   
/*ordena mejores de mayor a menor       
     
    
 */
  j=ordenarmejores();

tprom=0;
tprom=promediarmejores();
printf("----PromedioTotalV1:-- %1.5f",tprom);
//getch();
diversidaddinamica=0;
for(i=0;i<diversidad;i++){
   if(mejores[i][0]>tprom) diversidaddinamica++;
   };



 /* do recombination */
     for(i = 0; i < mu; i+= 2)
     {  
          if (drand(1) <= individual_recombination_probability)
          {
               if (variable_swap_probability > 0)
               {
                    result = uniform_crossover
                      (get_individual(result_ids[i]),
                       get_individual(result_ids[i + 1]));
                    if (result != 0)
                         log_to_file(log_file, __FILE__, 
                                     __LINE__, "recombination failed!");
               }

               if (variable_recombination_probability > 0)
               {
                    result = sbx
                      (get_individual(result_ids[i]), 
                       get_individual(result_ids[i + 1]));
                    if (result != 0)
                         log_to_file(log_file, __FILE__, 
                                     __LINE__, "recombination failed!");
               }
          }
     } 
     
     /* do mutation */
     for(i = 0; i < mu; i++)
     {
          if (drand(1) <= individual_mutation_probability)
          { 
               if (variable_mutation_probability > 0)
               {
                    result = mutation(get_individual(result_ids[i]));
                    if(result != 0)
                      log_to_file(log_file, __FILE__, __LINE__,
                                  "mutation failed!");
               }
          }
     } 
     

/* do Decision */        
for(i = 0; i < mu; i++)
     {
   if (drand(1) <= hill||hill==0.000123)
          { 
                 
               result = decision_making(get_individual(result_ids[i]));
               
                    if(result != 0)
                      log_to_file(log_file, __FILE__, __LINE__,
                                  "decision_making failed!");
               
          };
};     
//-------------------------------------------------------------     
  
     for(i = 0; i < mu; i++)
     {
          int result;
          result = eval(get_individual(result_ids[i]));
          
     }
    
     
     return (0);
}


int decision_making(individual *ind)
{ 
     int i,j,k,l,n1,k1,r;
     double r1,maxmutr1,retemp,reind,g,probind;
     individual *temp1=new_individual();
     reind=0;
     if (ind == NULL)
     {
	 return (1);
     }
     if(objetivo==1) g=1-(ind->f[0]);
            else g=1000-(ind->f[0]);
     if(ind->f[1]>0) reind=(g-p->rf)/sqrt(ind->f[1]);       
     for(j=0;j<p->numfondos;j++){
                 temp1->x[j]=ind->x[j];};
     i=0;
     k=0;
      l=0;
    
//----------------------------

if(reind<tprom&&i<diversidad){
k=irand(diversidaddinamica);
if(hill!=0.000123||(gen==0&&drand(1)<=hill2)){
for(j=0;j<p->numfondos;j++){
                 ind->x[j]=mejores[k][j+1];
                 temp1->x[j]=mejores[k][j+1];
                 //temp1->x[j]=maxxsharpe[j];
                 
                 };
                 reind=mejores[k][0];
                 retemp=mejores[k][0];
                 printf("%s","CAMBIADO POR ESTAR POR DEBAJO DEL PROMEDIO");
                 //getch();

l=diversidad;
k=0;
}
}else{insertarmejores(ind,reind);};

tprom=0;
tprom=promediarmejores();
diversidaddinamica=0;
for(i=0;i<diversidad;i++){
   if(mejores[i][0]>tprom) diversidaddinamica++;
   
   };
i=0;               
//printf("----diversidad dinamicaV1:-- %1i",diversidaddinamica);
if(gen==0&&drand(1)<=hill2){
probind=hill3;
printf("----CORRECTOR DE POBLACION INICIAL");
//getch();
}else{
probind=hill1;	
//printf("probind: %1.5f",probind);
//printf("probind: %1i",l);
}
if((drand(1)<=probind&&reind>tprom)||(l>0&&((gen>iteracion2*0.9)&&(gen<iteracion2)))){  
     j=irand(p->numfondos-1);
     k=0;
     //maxmutr1=mutation1(temp1,j);
     maxmutr1=drand(0.09);
    while(i<iteracion1){           
          r1 = temp1->x[j];
          nrocorridasBLV=0;
          if(r1>p->minfondo&&maxmutr1<p->maxfondo){ 
                temp1->x[j]=maxmutr1;
                retemp=0;
                nrocorridasBLV=1;
	              eval(temp1);
	              nrocorridasBL++;
	              //printf("----Objetivo:-- %1.5f",temp1->f[0]);
	              //getch();
	              //&&retemp>tprom &&reind>tprom
                  if(objetivo==1) g=1-(temp1->f[0]);
                   else g=1000-(temp1->f[0]);
                   if(temp1->f[1]>0) retemp=(g-p->rf)/sqrt(temp1->f[1]);
                   if(retemp>reind){
                      for(n1=0;n1<p->numfondos;n1++){               
                      ind->x[n1]=temp1->x[n1];};
                      reind=retemp;
                      maxmutr1+=drand(0.009);
                      printf("%s","----CAMBIADO POR HILL----");
                       insertarmejores(temp1,retemp);
                       tprom=0;
                      tprom=promediarmejores();
                       i=i+1;
                       k=1;
                      //getch();
                     }else{
                      if(reind>0){
                      if(((reind-retemp)/reind)>0.02){
                                 //printf("%s","NO.......CAMBIADO POR HILL");
                                 //maxmutr1=mutation1(temp1,j);
                                 maxmutr1=drand(0.9);
                                 i=i+iteracion1*0.3;
                                 //getch();
                      }else{
                          //maxmutr1=drand(p->maxfondo-(maxmutr1+(0.09)));
                          //maxmutr1=mutation1(temp1,j);
                          maxmutr1=drand(0.09);
                          i=i+iteracion1*0.1;
                          //printf("PASO AQUI");
                          };
                      i=i+iteracion1*0.1; 
                      if(i>=iteracion1) k=1;
                      for(n1=0;n1<p->numfondos;n1++){               
                      temp1->x[n1]=ind->x[n1];};
                      retemp=reind;
                      //printf("----ReindNOCAMBIADOTotal:-- %1.5f",retemp);
                      };
                 };//endif
          }else{
                  //printf("QUEDO AQUI");
                  i=i+iteracion1*0.1;
                  j=irand(p->numfondos);
           };//endif
    };//endwhile
  };//endif
    nrocorridasBLV=0;
    return (0);
}



int mutation(individual *ind)
{
     int i;

     if (ind == NULL)
     {
	 return (1);
     }
     
     for (i = 0; i < ind->n; i++)
     {
	 if (drand(1) <= variable_mutation_probability)
	 {
	     double eta = eta_mutation;
	     double u = drand(1.0);
	     double delta = 0;
	     double x = ind->x[i];
	     double lb = 0;    /* lower bound of variable i */
	     double ub = 1;    /* upper bound of variable i */
	     double diff = ub - lb;  /* range of variable i */
	     double maxmut0 = x - lb;
	     double maxmut = ub - x;
	     double delta_max = maxmut0 / diff;
	     if (maxmut0 > maxmut)
	     {
		 delta_max = maxmut / diff;
	     }
	     
	     if (u < 0.5)
	     {
		 double b =  2*u + (1-2*u)*(pow(1-delta_max,(eta+1)));
		 delta = pow(b,(1.0/(eta+1))) - 1.0;
	     }
	     else
	     {
		 double b = 2*(1-u) + 2*(u-0.5)*(pow(1-delta_max,(eta+1)));
		 delta = 1.0 - pow(b,(1.0/(eta+1)));
	     }
	     if (delta > delta_max)  /* machine accuracy problem */
		 delta = delta_max;
	     else if (delta < -delta_max)
		 delta = -delta_max;
	     
	     ind->x[i] = x + delta * diff;
	 }
     }
     
     return (0);
}

double mutation1(individual *ind,int i)
{
     //int i=irand(ind->n);

     if (ind == NULL)
     {
	 return (1);
     }
     
     if(ind->x[i]>0)
     {
	     double eta = eta_mutation;
	     double u = drand(1.0);
	     double delta = 0;
	     double x = ind->x[i];
	     double lb = 0;    /* lower bound of variable i */
	     double ub = 1;    /* upper bound of variable i */
	     double diff = ub - lb;  /* range of variable i */
	     double maxmut0 = x - lb;
	     double maxmut = ub - x;
	     double delta_max = maxmut0 / diff;
	     if (maxmut0 > maxmut)
	     {
		 delta_max = maxmut / diff;
	     }
	     
	     if (u < 0.5)
	     {
		 double b =  2*u + (1-2*u)*(pow(1-delta_max,(eta+1)));
		 delta = pow(b,(1.0/(eta+1))) - 1.0;
	     }
	     else
	     {
		 double b = 2*(1-u) + 2*(u-0.5)*(pow(1-delta_max,(eta+1)));
		 delta = 1.0 - pow(b,(1.0/(eta+1)));
	     }
	     if (delta > delta_max)  /* machine accuracy problem */
		 delta = delta_max;
	     else if (delta < -delta_max)
		 delta = -delta_max;
	     
	     ind->x[i] = x + delta * diff;
	 
     }
     
     return (ind->x[i]);
}

int uniform_crossover(individual *ind1, individual *ind2)
{
     int i;
  
     for (i = 0; i < ind2->n; i++)
     {
	 if (drand(1) <= variable_swap_probability) /* switch variable */
	 {
	     double x = ind2->x[i];
	     ind2->x[i] = ind1->x[i];
	     ind1->x[i] = x;
          } 
     }  

     return (0);
}

int unif_cross_memetic(individual *ind1, individual *ind2)
{
     int i;
  
     for (i = 0; i < ind2->n; i++)
     {
	 if (drand(1) <= variable_swap_probability) /* switch variable */
	 {
	     double x = ind2->x[i];
	     ind2->x[i] = ind1->x[i];
	     ind1->x[i] = x;
          } 
     }  

     return (0);
}

int sbx(individual *ind1, individual *ind2)
{
     int i;
  
     for (i = 0; i < ind2->n; i++)
     {
	 if (drand(1) <= variable_recombination_probability)  
	 {
	     double di = eta_recombination; /* distribution index */
	     int bounded = 1;
	     double lb = 0;    /* lower bound of variable i */
	     double ub = 1;    /* upper bound of variable i */	     
	     double u = drand(1);
	     double b0=0, b1=0;   /* spread factors */
	     double x0 = ind1->x[i];
	     double x1 = ind2->x[i];
	     double bl=0, bu=0, p_bl=0, p_bu=0, bll=0, buu=0, blll=0, buuu=0;
	     double dx = 0;
	     double u0=0, u1=0;
             
            /* calculate spread factor(s) b0, b1 */ 
            if (bounded == 1)
            {
                dx = fabs(x1-x0);   /* difference of x values */
                if (dx > 0)
                {
                    bl = (x0 + x1 - 2*lb) / dx;
                    bu = (2*ub - x0 - x1) / dx;
                    bll = (x0 + x1 - 2*(x0-lb)) / dx;
                    buu = (2*(ub-x1)-x0-x1) / dx;
                    if (x0 < x1)
                    {
                        blll = 1 + 2 * (x0 - lb) / dx;
                        buuu = 1 + 2 * (ub - x1) / dx;
                    }
                    else
                    {
                        blll = 1 + 2 * (x1 - lb) / dx;
                        buuu = 1 + 2 * (ub-x0) / dx;
		    }
		    
		    bl = blll; /* take Deb's version (numerically stable) */
                    bu = buuu;

                    /* switch off symmetric recombination to avoid
                     * getting stuck on a line where one parameter
                     * equals an extreme value. */
                    if (use_symmetric_recombination)
                    {
                       if (bl < bu)  /* symmetric distribution, like Deb */
                          bu = bl;
                       else
                          bl = bu;
                       assert (b0 == b1);
                    }
                    assert(bl > 0 && bu > 0);
                    p_bl = 1 - 1/(2*pow(bl,di+1));
                    p_bu = 1 - 1/(2*pow(bu,di+1));
                }
                else
                {
                    p_bl = 1;
                    p_bu = 1;
                }
                u0 = u*p_bl;
                u1 = u*p_bu;
                if (u0<=0.5)
                    b0 = pow(2*u0,1/(di+1));
                else
                    b0 = pow(0.5/(1-u0),1/(di+1));
                if (u1<=0.5)
                    b1 = pow(2*u1,1/(di+1));
                else
                    b1 = pow(0.5/(1-u1),1/(di+1));
                assert(dx==0 || (b0<=bl && b1<=bu)); /* machine accuracy */
            }
            else
            {
                if (u<=0.5)
                    b0 = pow(2*u,1/(di+1));
                else
                    b0 = pow(0.5/(1-u),1/(di+1));
                b1 = b0;
            }

            if (x0 < x1)
            {
                ind1->x[i] = 0.5*(x0+x1 + b0*(x0-x1));
                ind2->x[i] = 0.5*(x0+x1 + b1*(x1-x0));
            }
            else
            {
                ind1->x[i] = 0.5*(x0+x1 + b1*(x0-x1));
                ind2->x[i] = 0.5*(x0+x1 + b0*(x1-x0));
            }
	 }
     }  
     
     return (0);
}



/* Generate a random integer. */
int irand(int range)
{
     int j;
     j=(int) ((double) range * (double) rand() / (RAND_MAX + 1.0));
     return (j);
}
/*RETORNA EL COCIENTE DE UNA DIVISION*/
int cociente(int x,int y){
    div_t d;
    d=div(x,y);
return(d.quot);
};
    
/* Generate a random double. */
double drand(double range)
{
     
     double j;
     j=(range * (double) rand() / (RAND_MAX + 1.0));
     return (j);
}


void Leerfondos(ProblemInfo *p,char problema[FILE_NAME_LENGTH]){
     FILE *fp;
     int j,k,i;
     double r,aux,suma;
     char str[CFG_NAME_LENGTH];
     
       
     
     fp = fopen(problema, "r"); 
     assert(fp != NULL);
     
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "problem") == 0);
     fscanf(fp, "%s", p->problem);
     printf("%s",p->problem);
     fscanf(fp, "%s", str);
     assert(strcmp(str, "NumFondos") == 0);
     fscanf(fp, "%d", &(p->numfondos));
     printf("%d",p->numfondos);
      fscanf(fp, "%s", str);
     assert(strcmp(str, "NumBits") == 0);
     fscanf(fp, "%d", &(p->numbits));
      fscanf(fp, "%s", str);
     assert(strcmp(str, "varianza") == 0);
//     printf("%s",str);
           
      p->beneficios=(double *) malloc(sizeof(double)*p->numfondos);
     //busqueda local 
      p->sharpe=(double **) malloc(sizeof(double)*p->numfondos);
      for (j=0; j<p->numfondos; j++)
          p->sharpe[j]=(double *) malloc(sizeof(double)*p->numfondos);
     //busqueda local
     p->varfondos=(double **) malloc(sizeof(double)*p->numfondos);
      for (j=0; j<p->numfondos; j++)
          p->varfondos[j]=(double *) malloc(sizeof(double)*p->numfondos);
     //mejores para los variadores
      /*mejores=(double **) malloc(sizeof(double)*(p->numfondos+1));
      for (j=0; j<p->numfondos+1; i++)
          mejores[j]=(double *) malloc(sizeof(double)*p->numfondos);*/
    
      for (j=0; j<p->numfondos; j++){
        for (k=0; k<p->numfondos; k++){
           fscanf(fp, "%lf",&r);
           p->varfondos[j][k]=r;

           };
            };
      
      fscanf(fp, "%s", str);
     assert(strcmp(str, "beneficio") == 0);
     for (j=0; j<p->numfondos; j++){
       fscanf(fp, "%lf",&r);
       p->beneficios[j]=r;};
     fscanf(fp, "%s", str);
     assert(strcmp(str, "Rf") == 0);
     fscanf(fp, "%lf", &r);
     p->rf=r;
  
    fscanf(fp, "%s", str);
    assert(strcmp(str, "Rfmax") == 0);
     fscanf(fp, "%lf", &r);
     p->rfmax=r;
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "Rfmin") == 0);
     fscanf(fp, "%lf", &r);
     p->rfmin=r;
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "Cardinalidad") == 0);
     fscanf(fp, "%d", &(p->cardinalidad));
    
     fscanf(fp, "%s", str);
     assert(strcmp(str, "maxfondo") == 0);
     fscanf(fp, "%lf", &r);
     p->maxfondo=r;
     
     fscanf(fp, "%s", str);
     assert(strcmp(str, "minfondo") == 0);
     fscanf(fp, "%lf", &r);
     p->minfondo=r;
     printf("%s","Termino de leer");
     fclose(fp);
     
    

}


/* Determines the objective value based on DTLZ */
double eval(individual *ind)
{
if(nrocorridasBLV==0){
   nrocorridas=nrocorridas+1.0;
}else{
      nrocorridas=nrocorridas+(2.0 / ((double)p->numfondos+1.0));
      //printf("----Corrida:-- %f ---%d----%d",nrocorridas,nrocorridasBLV,nrocorridasBL);
      //getch();
      };
   if (strcmp(p->problem, "VARIABLE") == 0)
    {
	return (Evaluar(ind,p));
    }
    if (strcmp(p->problem, "MIXTO") == 0)
    {
    return (Evaluar(ind,p));
    }
    if (strcmp(p->problem, "FIJO") == 0)
    {
	return (Evaluar(ind,p));
    }
     if (strcmp(p->problem, "EUROPA") == 0)
    {
	return (Evaluar(ind,p));
    }
      log_to_file(log_file, __FILE__, __LINE__, "unknown problem specified");
    return (1);
}


/*--------------------------GASSHARPE---------------------------*/
void AnalizaRestricciones (double* pesos, ProblemInfo* p2)
{
       int indice[100];
       int i,j;
       double value, total=0;
       
       for (i=0; i<p2->numfondos; i++)
           indice[i] = i+2;

/*       for (i=0; i<p->numfondos; i++)
           printf ("%f[%d] ", pesos[indice[i]], indice[i]);
       printf ("\n"); 
*/
       for (i = 0; i<p2->numfondos; i++) 
            for (j = i+1; j<p2->numfondos; j++) 
               if (pesos[indice[j]]>pesos[indice[i]]) {
                  int tmp = indice[i];
                  indice[i] = indice[j];
                  indice[j] = tmp;
               }
/*       for (i=0; i<p->numfondos; i++)
           printf ("%f[%d] ", pesos[indice[i]], indice[i]);
       printf ("\n"); 
*/ 
 
           
       for (i=0; i<p2->cardinalidad; i++) {
           if (pesos[indice[i]]==0)
              pesos[indice[i]]=1;
           total += pesos[indice[i]];
       }
       for (i=p2->cardinalidad; i<p2->numfondos; i++)
           pesos[indice[i]] = 0;
/*       printf ("Prenormalizado:");
       for (i=0; i<p->numfondos; i++)
           printf ("%f[%d] ", pesos[indice[i]], indice[i]);
       printf ("\n"); */
  //if ((p2->cardinalidad*p2->minfondo<1)&&(p2->cardinalidad*p2->minfondo>0)){
     value = (1.0-p2->cardinalidad*p2->minfondo)/total;
       for (i=0; i<p2->cardinalidad; i++){
           pesos[indice[i]] = p2->minfondo + pesos[indice[i]]*value;
          
              };
//};   
 
 
 
 
}

/*
 * Evalua un individuo.
 *
 */
 
double Evaluar (individual *ind, void* info){
    int j,k,r,indice1[100];
    double total, beneficio, riesgo, f1, f2, total1,total2,total3,ft;
   ProblemInfo *p1= (ProblemInfo *) info;
    double* fenotipo;
     fenotipo = (double*) malloc (sizeof(double)*(p1->numfondos+2));
     
     
    r=0;
    total=0;
        
    for (j=0; j<p1->numfondos; j++){
        fenotipo[j+2]=0;
        fenotipo[j+2]=ind->x[r];
        r++;
        total+= fenotipo[j+2];
          };
        
  //printf("Cardinalidad: %lf",p1->cardinalidad);
if (p1->cardinalidad>0){
   
       AnalizaRestricciones (fenotipo, p1);
       
/*---------------- MAXFONDO------------------------------------     */     
              total1=0;
              total=0;
              total2=0;
              total3=0;
          if (p1->maxfondo>0){
             if (p1->cardinalidad*p1->maxfondo>=1){
                 for (j=0; j<p1->numfondos; j++){
                     if (fenotipo[j+2]>p1->maxfondo){
                          fenotipo[j+2]=p1->maxfondo;
                          indice1[j]=1.00;
                          }
                     else indice1[j]=0.00;
                   total+=fenotipo[j+2];
                  }; 
                 if (total>=1){
                     for (j=0; j<p1->numfondos; j++){  
                               fenotipo[j+2]/=total;
                                              };
                      }
                 else {
                     for (j=0; j<p1->numfondos; j++){
                          if (indice1[j]==1) total1+= fenotipo[j+2];
                          else total2+= fenotipo[j+2];
                      };
                      total=0;
                      total3=0;
                      for (j=0; j<p1->numfondos; j++){
                          if (indice1[j]==0){ 
                             fenotipo[j+2]+=fenotipo[j+2]*(1.000-total1)/total2;
                             if (fenotipo[j+2]>p1->maxfondo) fenotipo[j+2]=p1->maxfondo;
                             else total3+=fenotipo[j+2];
                             };
                          total+=fenotipo[j+2];      
                       /*for (j=0; j<p->numfondos; j++)
                           printf ("%f", fenotipo[j+2]);
                       getchar();*/
                     
                      };
                      for (j=0; j<p->numfondos; j++){
                           if (total>1) fenotipo[j+2]/=total;
                           else {
                                if (fenotipo[j+2]<p1->maxfondo) fenotipo[j+2]+=fenotipo[j+2]*(1.000-total)/total3;};};
                 };
             };
          };        
    /*---------------- MAXFONDO------------------------------------*/  
         r=0;
        total=0;
        for (j=0; j<p1->numfondos; j++){
          ind->x[r]=fenotipo[j+2];
          r++;
          total+= fenotipo[j+2];
          };
        /*for (j=0; j<p1->numfondos; j++){
              ind->x[j]/=total;
              fenotipo[j+2]=ind->x[j];}; */  
           
       }
    else
        for (j=0; j<p1->numfondos; j++){
            ind->x[j]/=total;
            fenotipo[j+2]=ind->x[j];
            };
 
 
 
           
            
    fenotipo[0]=0;
    for (j=0; j<p1->numfondos; j++)
        fenotipo[0]+=fenotipo[j+2]*p1->beneficios[j];
    
    fenotipo[1]=0;
    for (j=0; j<p1->numfondos; j++){
        for (k=0; k<p1->numfondos; k++){
            fenotipo[1]+=fenotipo[j+2]*fenotipo[k+2]*p1->varfondos[j][k];};
        }
    
    
if (objetivo==1){            
    f1=-1*fenotipo[0];   
    
    if (1+f1>0) ind->f[0] = 1+f1;
    else ind->f[0] = 999;  
    
    
    if  (fenotipo[1]>0) ind->f[1] = fenotipo[1];
    else ind->f[1] = 999.99; 
    ft=(ind->f[0]-p->rf)/sqrt(ind->f[1]);  
}else{
    
    f2= -1*(fenotipo[0]-p->rf)/sqrt(fenotipo[1]);
    if (1000+f2>0) ind->f[0] = 1000+f2;
    else ind->f[0] = 999;  
    
    if  (fenotipo[1]>0) ind->f[1] = fenotipo[1];
    else ind->f[1] = 999.99; 
    ft=ind->f[0];
};
  
  
    return(ft);
}
/*-----------------------HASTA AQUI!!!------------------------------*/



/* create a random new individual and allocate memory for it,
   returns a pointer to the new individual */
individual *new_individual()
{
     individual *return_ind;
     int i;
     int k=Card;
     int g=0;
    
     

     return_ind = (individual *) malloc(sizeof(individual));
     return_ind->x = (double *) malloc(sizeof(double) * p->numfondos);
     return_ind->f = (double *) malloc(sizeof(double) * dimension);
     return_ind->h = (double *) malloc(sizeof(double) * p->numfondos);
     /*-----------------GASSHARPE-----------------------*/
     return_ind->Phenotype = (double*) malloc (sizeof(double)*(p->numfondos+2));
     for (i = 0; i < p->numfondos; i++)
     {
	 return_ind->h[i] = 0.0;
     }
     i=0;
     do{
      i=irand(p->numfondos);
      if(return_ind->h[i]==0){
            return_ind->h[i]=1.0;
            g++;
            };
      }while(g<k);
      
     for (i = 0; i < p->numfondos; i++)
     {
	 //return_ind->x[i] =return_ind->h[i]*drand(randcard);
	 
	 return_ind->x[i] =drand(randcard);
	 return_ind->Phenotype[i+2]=return_ind->x[i];
     }
     
     return_ind->n = p->numfondos;

     for (i = 0; i < dimension; i++)
     {
	 return_ind->f[i] = 0.0;
	 return_ind->Phenotype[i]=0.0;
     }
     
     return (return_ind);
}


/* copy an individual and return the pointer to it */
individual *copy_individual(individual *ind)
{
     individual *return_ind;
     int i;
     
     return_ind = (individual *) malloc(sizeof(individual));
     return_ind->x = (double *)malloc(sizeof(double) * p->numfondos);
     return_ind->f = (double *) malloc(sizeof(double) * dimension);
     return_ind->h = (double *)malloc(sizeof(double) * p->numfondos); 
/*-----------------GASSHARPE-----------------------*/
     return_ind->Phenotype = (double*) malloc (sizeof(double)*(p->numfondos+2));     
     for (i = 0; i < p->numfondos; i++)
          return_ind->x[i] = ind->x[i];

     for (i = 0; i < dimension; i++)
	  return_ind->f[i] = ind->f[i];

     return_ind->n = ind->n;
     
for (i = 0; i < p->numfondos; i++)
          return_ind->h[i] = ind->h[i];
    /*-----------------GASSHARPE-----------------------*/
for (i = 0; i < p->numfondos+dimension; i++)
          return_ind->Phenotype[i] = ind->Phenotype[i];  
     return(return_ind);
}


int log_to_corrida(char *file)
{
     FILE *fp;
       
     if(file != NULL)
     {
          fp= fopen(file, "a");
          assert(fp != NULL);
          if(hill>0){ fprintf(fp, "%sH  %i   %d  %d", filenamebase_selector,nrocorridasmax,nrocorridasBL,impgen);
          }else fprintf(fp, "%s  %f", filenamebase_selector,nrocorridas);
          fprintf(fp, "\n");
          nrocorridas=0.0;
          nrocorridasBL=0;
          nrocorridasBLV=0;
          };
          
            

          fclose(fp);
    
     return (0);
}
//Vector mejores Sharpe---------------------------
int ordenarmejores()
{
int j,k,nfondo,n;
double burbuja[50];

 for (j=0; j<diversidad; j++){
           for(k=j+1;k<diversidad;k++){                       
              if (mejores[j][0]<mejores[k][0]){
                   for(n=0;n<50;n++) burbuja[n]=0;
                   burbuja[0]=mejores[j][0];
                   for(nfondo=1;nfondo<p->numfondos+1;nfondo++){
                      burbuja[nfondo]=mejores[j][nfondo];
                   };
                   mejores[j][0]=mejores[k][0];
                   for(nfondo=1;nfondo<p->numfondos+1;nfondo++){
                      mejores[j][nfondo]=mejores[k][nfondo];
                   };
                   mejores[k][0]=burbuja[0];
                  for(nfondo=1;nfondo<p->numfondos+1;nfondo++){
                      mejores[k][nfondo]=burbuja[nfondo];
                   };
              };
           };         
       };    
    return(0);
}


double promediarmejores()
{
       int i;
       double tprom1=0.0;
  for(i=0;i<diversidad;i++){
   //printf("----MEJORES SharpePromedio1:-- %1.10f",mejores[i][0]); 
   tprom1+=mejores[i][0];
   
   };
  //getch();
tprom1=tprom1/diversidad;     
return(tprom1);      
       }
       
int insertarmejores(individual *ind, double sharpe){
 int ver;
 int j=0;
 int k=0;
 int k1=0;
 int r12,nfondo;
             
               for(r12=0;r12<diversidad;r12++){
                        if(mejores[r12][0]==sharpe) k1=1;                     
                                             };
               for(j=0;j<diversidad;j++){
                   //printf("mejores1: %1.5f",mejores[j][0]);
                  if(mejores[j][0]==0&&k1==0){ 
                   mejores[j][0]=sharpe;
                   //mejores[j][30]=g1;
                  printf("%s","LO METIO EN LOS CEROS");
                   //printf("J***: %1i",j);
                   for(nfondo=1;nfondo<p->numfondos+1;nfondo++){
                      mejores[j][nfondo]=ind->x[nfondo-1];
                   };
                   break;
                  };      
              };
              k1=0;
              j=ordenarmejores();
                       
               for(r12=0;r12<diversidad;r12++){
                        if(mejores[r12][0]==sharpe) k1=1;                     
                                             };            
               for(j=0;j<diversidad;j++){
                if(mejores[j][0]<sharpe&&k1==0){
                  printf("----SharpeM:-- %1.10f",sharpe);     
                  //getch();      
                   mejores[diversidad-1][0]=sharpe;
                   for(nfondo=1;nfondo<p->numfondos+1;nfondo++){
                      mejores[diversidad-1][nfondo]=ind->x[nfondo-1];
                   };
                   ordenarmejores();
                     /*for(ver=0;ver<diversidad;ver++){
                    printf("----MEJORES SharpeDespues:-- %1.10f",mejores[ver][0]); };
                    getch();*/
                      break;
                  };
               };
 return(k);
           }
           
           
/* Writes the index, objective values and bit string of
   all individuals in global_population to 'out_filename'. */
void write_output_file()
{
     int j, current_id;
     /*int i = 0;
     int n = p->NumFondos;*/
     double g = 0;
     FILE *fp_out;
     individual *temp;
     char outfile1[FILE_NAME_LENGTH];
        
     if(hill>0||hill2>0){
     sprintf(outfile1, "%sH%s%sE%d%K%dPM%iP%i.%d", filenamebase_selector,p->problem,modelo,currentRun+1,p->cardinalidad,(int)(hill*100),(int)(hill1*100),maxgen);
     }else sprintf(outfile1, "%s%s%sE%d%K%d.%d", filenamebase_selector,p->problem,modelo,currentRun+1,p->cardinalidad,maxgen);
     
     
     
     
     fp_out = fopen(outfile1, "w");
     assert(fp_out != NULL);
     current_id = get_first();

     while (current_id != -1)
     {       
	  temp = get_individual(current_id);
	  
	 
      
          fprintf(fp_out, "%d ", current_id); /* write index */
	  for (j = 0; j < dimension; j++)
	  {
	      if (j==0) {
            if(objetivo==1) g=1-(temp->f[j]);
            else g=1000-(temp->f[j]);
            fprintf(fp_out, "%f ", g);}
          else {
               /*if (j==2) {
               g=1000-(temp->f[j]);
               fprintf(fp_out, "%f ", g);}
               else {*/fprintf(fp_out, "%f ", temp->f[j]);
              //};
               };
	  }
          
	  for (j = 0; j < temp->n; j++)
          {
               printf("%f",temp->x[j]);
               fprintf(fp_out, "%f ", temp->x[j]);
          }
          fprintf(fp_out, "\n");
   /*}*/
	  current_id = get_next(current_id);
     }
     printf("%s","llego para copiar archivo");
     updateVariatorSeed();
     printf("%s","Termino la semilla");
     fclose(fp_out);
     
     log_to_corrida(nro_corridas);
     
     
     printf("%s",paramfile);
}

void updateVariatorSeed(){

	FILE *fp = NULL;
	char *(contents[LINE_LENGTH]);
	int linectr = 0;
	int line;
	char str[CFG_NAME_LENGTH]; 
	int oldSeed, newSeed;

	printf("    update variator seed\n");

	// read from parameter file into array of strings contents
printf("%s",paramfile);
	fp = fopen(paramfile, "r");
	assert(fp != NULL);

	contents[0] = (char *) malloc(sizeof(char[LINE_LENGTH]));
	assert(contents[0] != NULL);

	while (fgets(contents[linectr], LINE_LENGTH, fp) != NULL) {
		linectr++;
		contents[linectr] = (char *) malloc(sizeof(char[LINE_LENGTH]));
		assert(contents[linectr] != NULL);
	}
	free(contents[linectr]);
	fclose(fp);

	// write to parameter file from array of strings contents
	// replace seed

	fp = fopen(paramfile, "w");
	assert(fp != NULL);

	line = 0;
	while (line < linectr) {
		sscanf(contents[line], "%s %d", str, &oldSeed);
		if (strcmp(str, "seed") == 0) {
			newSeed = irand(RAND_MAX);
			fprintf(fp, "seed %d\n", newSeed);
		printf("      seed %d replace by %d\n", oldSeed, newSeed);
		} else {
			fputs(contents[line], fp);
		}
		free(contents[line]);
		line++;
	}
	fclose(fp);
}



/**********| addition for DTLZ end |*******/
