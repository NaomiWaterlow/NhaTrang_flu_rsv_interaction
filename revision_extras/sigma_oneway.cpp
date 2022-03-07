#include <Rcpp.h>
#include <R.h>
using namespace Rcpp;

//This is your global parms variable that can be accessed in both initmod and derivs (or other functions)
Rcpp::List global_parms;

NumericVector transmission_calc_R2(double *y,
                                   NumericMatrix Contact_Structure,
                                   int num_grps, float bR){
  //  Rprintf("\n y0 is %f ", y[0]);
  //initialise output vector
  NumericVector transmission_R(num_grps*14);
  //for the age group being infected
  for(int l=0;(num_grps)>l;++l){
    //set up as 0 initially
    transmission_R[l]=0;
    //for each infecting age group
    for (int k=0;(num_grps)>k;++k){
      //wotk out relative infection
      float temp = bR * ( y[1+k*14] * Contact_Structure(l, k) + y[5+k*14] * Contact_Structure(l, k) +
                          y[8+k*14] * Contact_Structure(l, k));
      //add to previous infections of age group
      transmission_R[l]= transmission_R[l]+temp;
    };
    //   Rprintf("\n y1 is %f ", y[1]);
  }; return(transmission_R);
}

NumericVector transmission_calc_I2 (double *y,
                                    NumericMatrix Contact_Structure,
                                    int num_grps, float bI){
  //initialise output vector
  NumericVector transmission_I(num_grps*14);
  //for the age group being infected
  for(int a=0;(num_grps)>a;++a){
    //set up as 0 initially
    //  Rprintf("%i a", a);
    transmission_I[a]=0;
    //Rprintf("\n here 2");
    //for each infecting age group
    for (int b=0;(num_grps)>b;++b){
      //wotk out relative infection
      float temp = bI * ( y[4+b*14] * Contact_Structure(a, b) + y[5+b*14] * Contact_Structure(a, b) +
                          y[6+b*14] * Contact_Structure(a, b));
      //add to previous infections of age group
      transmission_I[a]= transmission_I[a]+temp;
      // Rprintf("\n transmission I is %f and temp is %f",transmission_I[a],temp);
    };
  }; return(transmission_I);
}

//Function to return the parms list passed to deSolve as a SEXP object
SEXP attribute_hidden get_deSolve_gparms_Rcpp()
{
  static SEXP(*fun)() = NULL;
  if (fun == NULL)
    fun = (SEXP(*)()) R_GetCCallable("deSolve","get_deSolve_gparms");
  return fun();
}
// 
extern "C" {
  void derivatives (int *neq, double *t, double *y, double *ydot,
                    double *yout, int *ip);
  void initmod(void(* odeparms)(int *, double *));
}

void initmod (void(* odeparms)(int *, double *))
{
  //We get the parms argument passed to deSolve as SEXP object
  SEXP sparms = get_deSolve_gparms_Rcpp();
  try {
    //Parse parameters passed to deSolve as Rcpp::List
    Rcpp::List parms = Rcpp::clone(Rcpp::as<Rcpp::List>(sparms));
    //We can't do garbage collection at end of model run in deSolve, so we do it if the same DLL/SO is still loaded
    // and deSolve is ran again
    global_parms = parms;
  } catch(std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
}

void derivatives (int *neq, double *t, double *y, double *ydot,
                  double *yout, int *ip)
{
  
  int num_grps = global_parms["num_grps"];
  float bR = global_parms["betaR"];
  float bI = global_parms["betaI"];
  float sigma_R = global_parms["sig_fixed"];
  float sigma_I = global_parms["sig"];
  float gammaR =  global_parms["gammaR"];
  float gammaI = global_parms["gammaI"];
  float rho =  global_parms["rho"];
  double seedI = global_parms["seedI"];
  double seedR = global_parms["seedR"];
  float trickleI_t = global_parms["trickleI"];
  float trickleR_t = global_parms["trickleR"];
  NumericVector RSV_Sus = global_parms["RSV_Sus"];
  NumericMatrix Contact_Structure = global_parms["Contact_Structure"];
  NumericVector trickleI(5);
  NumericVector trickleR(5);
  
  
  NumericVector N = global_parms["pop_sizes"];
  //  Rprintf("how many time this shown");
  if (*t<seedI){ 
    trickleI=rep(0,5);
  } else {
    trickleI=trickleI_t*N;}
  
  if (*t<seedR){
    trickleR=rep(0,5);
  } else {
    trickleR=trickleR_t*N;}
  
  NumericVector transmission_R = transmission_calc_R2 (y, Contact_Structure,num_grps,bR);
  NumericVector transmission_I = transmission_calc_I2 (y, Contact_Structure,num_grps,bI);
  
  for( int i=0; (num_grps)>i; ++i){
    
    //  Rprintf("\n time  %f ",t);
    
    //SS
    ydot[0+i*14]= ( - (RSV_Sus[i]*transmission_R[i]*y[0+i*14])
                      - (transmission_I[i]*y[0+i*14])
                      -trickleI[i]
                      -trickleR[i]
    );
    
    
    //     Rprintf("\n at time %f trickleI is %f", t, trickleI);
    
    
    
    //IS
    ydot[1+i*14]= ( (RSV_Sus[i]*transmission_R[i]*y[0+i*14])
                      - (transmission_I[i]*sigma_I*y[1+i*14])
                      - (gammaR*y[1+i*14])
                      + trickleR[i]
    );
    
    //PS
    ydot[2+i*14] = ( (gammaR*y[1+i*14])
                       - (rho*y[2+i*14])
                       - (sigma_I*transmission_I[i]*y[2+i*14]) );
                       
                       //RS
                       ydot[3+i*14] = ( (rho*y[2+i*14])
                                          - (transmission_I[i]*y[3+i*14]));
                       
                       
                       //SI
                       ydot[4+i*14] = ( (transmission_I[i]*y[0+i*14])
                                          - (RSV_Sus[i]*transmission_R[i]*sigma_R*y[4+i*14])
                                          - (gammaI *y[4+i*14])
                                          + trickleI[i]
                       )
                         ;
                       
                       //   Rprintf("\n at time %f transmissionR is %f", t, transmission_R[1]);
                       //   Rprintf("\n at time %f RSV_SUS[i] is %f", t,(RSV_Sus[0]) ); _FINE
                       //  Rprintf("\n at time %f gamaI is %f", t,gammaI );
                       // Rprintf("\n at time %f 1st line is %f", t, transmission_I[0]*y[0+0*14]);
                       
                       
                       
                       
                       //II
                       ydot[5+i*14] = ( (transmission_I[i]*y[1+i*14]*sigma_I)
                                          + (RSV_Sus[i]*transmission_R[i]*sigma_R*y[4+i*14])
                                          - (gammaR*y[5+i*14])
                                          - (gammaI*y[5+i*14]) ) ;
                                          //
                                          //P/RI
                                          ydot[6+i*14] = ( (transmission_I[i]*(sigma_I*y[2+i*14]+y[3+i*14]) )
                                                             + (gammaR*y[5+i*14])
                                                             - (gammaI*y[6+i*14]));
                                                             
                                                             
                                                             //SP
                                                             ydot[7+i*14] = ( (gammaI*y[4+i*14])
                                                                                - (rho*y[7+i*14])
                                                                                - (RSV_Sus[i]*sigma_R*transmission_R[i]*y[7+i*14]) );
                                                                                
                                                                                
                                                                                //IP/R
                                                                                ydot[8+i*14] = ((RSV_Sus[i]*transmission_R[i]*(sigma_R*y[7+i*14]+y[9+i*14]))
                                                                                                  + (gammaI*y[5+i*14])
                                                                                                  - (gammaR*y[8+i*14]));
                                                                                                  
                                                                                                  //SR
                                                                                                  ydot[9+i*14] = ( (rho*y[7+i*14])
                                                                                                                     - (RSV_Sus[i]*transmission_R[i]*y[9+i*14]) );
                                                                                                  
                                                                                                  //   Rprintf("\n at time %f transmissionR is %f", t, transmission_R[1]);
                                                                                                  //   Rprintf("\n at time %f RSV_SUS[i] is %f", t,(RSV_Sus[0]) ); _FINE
                                                                                                  //Rprintf("\n at time %f transmission_R[i] is %f", t,transmission_R[0] ); __FINE
                                                                                                  //     Rprintf("\n at time %f 7i is %f", t, - (rho*y[7+i*14]));
                                                                                                  
                                                                                                  //RR
                                                                                                  ydot[10+i*14] = ((gammaR*y[8+i*14])
                                                                                                                     + (gammaI*y[6+i*14]) );
                                                                                                  
                                                                                                  //R cases
                                                                                                  ydot[11+i*14] = ((RSV_Sus[i]*transmission_R[i]*y[9+i*14])
                                                                                                                     + (RSV_Sus[i]*sigma_R*transmission_R[i]*y[7+i*14])
                                                                                                                     + (RSV_Sus[i]*transmission_R[i]*y[0+i*14]) 
                                                                                                                     + trickleR[i]
                                                                                                  );
                                                                                                  
                                                                                                  //I cases
                                                                                                  ydot[12+i*14] =( (transmission_I[i]*y[3+i*14])
                                                                                                                     + (transmission_I[i]*sigma_I*y[2+i*14])
                                                                                                                     + (transmission_I[i]*y[0+i*14]))
                                                                                                                     + trickleI[i]
                                                                                                                     ;
                                                                                                  // Dual cases
                                                                                                  ydot[13+i*14] = ( (transmission_I[i]*y[1+i*14]*sigma_I)
                                                                                                                      + (RSV_Sus[i]*transmission_R[i]*sigma_R*y[4+i*14])
                                                                                                                      
                                                                                                  );
                                                                                                  
                                                                                                  // 
                                                                                                  // Rprintf("\n time  %f ",t);
                                                                                                  //          Rprintf("\n time  %f value  %.20f",t,(transmission_I[i]*y[1+i*14]*sigma_both));
                                                                                                  // 
                                                                                                  
  };
  // Rprintf("\n %f gives IS1 %f",t, ydot[1+0*14]);
  // return(ydot);
  // Rprintf("\n %f", y_out[12+5*14]);
  // ydot;
  // return List::create(_["yout"] = ydot) ;
  
}