#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
//#include <string>
//#include <assert.h>
#include <stdlib.h>
#include "generate_trajectories.h"
#include "derivs.h"
#include "rkf.h"
#include "essentials.h"
#include "prms.h"
#include "assert.h"


// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double* yic;
extern prms* ppc;
extern double G_CLO_INTRODUCTION_TIME;
extern int G_CLO_INTRODUCTION_COUNT;

//  END  ### ### GLOBAL VARIABLES ### ###

// [[Rcpp::export]]
void generate_trajectories( double inital_number_of_cases, double param_beta, double t0, double tf, double h )
{
    assert(ppc);
    double tt,rkstep,ttstop,eps,h1,hmin,ttbeforestop;
    int nvar,nok,nbad; //rkqs();
    int i, j;
    int dim = STARTJ+NUMAC; // this is the dimension of the dynamical system, i.e. the number of equations
    
    
    // NOTE 2020-04-04 : this is still very fast and it prevents an off-by-one-day pseudo-error when the step size is close to 0.5 or 1.0
    //                 : just to be clear, a step-size of 1.0 is perfectly fine, but you have to mentally correct for the fact that the 
    //                 : diff-eqs will be about one day late
    int steps_per_day = 100;
    
    // set some values in the RKF integrator
    tt = (double)((int)t0); // chop off any decimals, and just start at an integer time
    ttstop = (double)((int)tf); // chop off any decimals, and just start at an integer time
    rkstep = 1.0 / ((double)steps_per_day);
    eps=0.00000001;
    h1 = 0.1;
    hmin = 1e-13;
    //hmin =    0.00000000001;
    nvar = dim;
    
    int counter = 0;
    bool bIntroductionOccurred = false;
    double NextTimeForBetaUpdate=1000000.0; assert( tf < 999999.0 );
    if( ppc->v_betatimes.size() > 1 )
    {
        NextTimeForBetaUpdate = ppc->v_betatimes[1];
    }
    ppc->assign_new_beta(); // called here for the first time, so it just puts the first beta into use
    
    // before integration begins, assign the initial beta value
    
    //
    //BEGIN OF LOOP INTEGRATING ODEs
    //
    while( tt < ttstop )
    {
        // introduce the first infections
        if( !bIntroductionOccurred && tt > G_CLO_INTRODUCTION_TIME )
        {
            yic[STARTI+4]=((double)G_CLO_INTRODUCTION_COUNT); // some number of infected individuals in their forties are introduced 
            bIntroductionOccurred = true;
        }
        
        // check if the beta value needs to be updated
        if( tt > NextTimeForBetaUpdate )
        {
            ppc->assign_new_beta();
            NextTimeForBetaUpdate = ppc->get_new_update_time();
        }
        
        if ( odeint(yic,nvar,tt,tt+rkstep,eps,&h1,hmin,&nok,&nbad, derivs,rkqs) != 1)
      	{
	    fprintf(stderr, "\n\nEXITING BC OF ERROR IN ODEINT\n\n");
            exit(-1);
        }

        if( counter%steps_per_day==0 )
        {
            printf("\n%1.3f", tt);
            for(i=0;i<dim;i++) printf("\t%1.4f", yic[i]);
            //for(i=0;i<NUMAC+3;i++) printf("\t%1.4f", yic[i]);
        }
        
	tt += rkstep;

        counter++;
    }
    //
    //END OF LOOP INTEGRATING ODEs
    //

    //printf("\n\n    %d \n\n", Q );
    
    delete[] yic;    
    
    //return vv;
}

