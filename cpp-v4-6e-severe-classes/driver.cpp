#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <assert.h>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <fstream>
#include <vector>

#include "generate_trajectories.h"
#include "derivs.h"
#include "rkf.h"
#include "prms.h"
#include "essentials.h"
#include "parseargs.h"

#include <ctime>
#include <sys/time.h>

using namespace std;


// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double b,d,s,c0,c2,c3;	// this are defined, i believe, in rkf.cpp
double* yic;

FILE* OutFile = NULL;

prms* ppc;  // need to allocate space for this in main

double G_CLO_BETA = 0.3;
double G_CLO_INTRODUCTION_TIME;
int G_CLO_INTRODUCTION_COUNT;
string G_CLO_LOCATION = "RI";

double G_CLO_TF = 365.0;
double G_CLO_P_HOSP_TO_ICU = 0.30;

bool G_B_DIAGNOSTIC_MODE = false;

// END ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );


//
//
int main(int argc, char* argv[])
{
    

    //
    // ###  1.  ALLOCATE SPACE FOR A PARAMETERS CLASS AND FOR THE MAIN STATE VARIABLES
    //
    ppc = new prms; assert( ppc );
    yic = new double[STARTK+NUMAC];
    for(int i=0;i<STARTK+NUMAC;i++) yic[i]=0.0; // zero everything out


    // assign some default values to the v_beta and v_betatimes arrays
    ppc->v_betas.push_back( 1.0 );
    ppc->v_betas.push_back( 0.9 );
    ppc->v_betas.push_back( 0.8 );
    ppc->v_betatimes.push_back( 0.0 );
    ppc->v_betatimes.push_back( 60.0 );
    ppc->v_betatimes.push_back( 95.0 );
    
//     string s1(argv[0]); string s2(argv[1]); string s3("3.8.44"); 
//     
//     printf("\n\t %s float outcome %d", s1.c_str() , isFloat(s1)?1:0 );
//     printf("\n\t %s float outcome %d", s2.c_str() , isFloat(s2)?1:0 );
//     printf("\n\t %s float outcome %d", s3.c_str() , isFloat(s3)?1:0 );
//     printf("\n\n");
    
    
    ParseArgs( argc, argv );

    
    /*for(int i=0; i < ppc->v_betas.size(); i++)
    {
        printf("\n\t time-start: %1.3f   \t beta-val: %1.3f", ppc->v_betatimes[i], ppc->v_betas[i] );
    }
    printf("\n\t time final is %1.1f\n", G_CLO_TF);
    
   return 1;*/
    
    
    
    //
    // ###  2.  INITIALIZE PARAMETERS - these are the default/starting values
    //
    ppc->v[ i_len_incub_period ]                            = 6.0;
    ppc->v[ i_len_symptomatic_infectious_period_phase_1 ]   = 7.0;
    ppc->v[ i_len_symptomatic_infectious_period_phase_2 ]   = 7.0;
    ppc->v[ i_len_medicalfloor_hospital_stay ]              = 12.0;
    
    // set the fraction of individuals who are hospitalized immediately after I_2
    ppc->v_fraction_hosp[0] = 0.0001;
    ppc->v_fraction_hosp[1] = 0.0002;
    ppc->v_fraction_hosp[2] = 0.0002;
    ppc->v_fraction_hosp[3] = 0.0002;
    ppc->v_fraction_hosp[4] = 0.005;
    ppc->v_fraction_hosp[5] = 0.01;
    ppc->v_fraction_hosp[6] = 0.02;
    ppc->v_fraction_hosp[7] = 0.07;
    ppc->v_fraction_hosp[8] = 0.12;

    // set the fraction of individuals who are immediately admitted to critical care after I_2
    ppc->v_fraction_crit[0] = 0.0001;
    ppc->v_fraction_crit[1] = 0.0002;
    ppc->v_fraction_crit[2] = 0.0002;
    ppc->v_fraction_crit[3] = 0.0002;
    ppc->v_fraction_crit[4] = 0.0007;
    ppc->v_fraction_crit[5] = 0.01;
    ppc->v_fraction_crit[6] = 0.02;
    ppc->v_fraction_crit[7] = 0.07;
    ppc->v_fraction_crit[8] = 0.12;

    // set the probability of death for I4 for all 9 age classes
    ppc->v_prob_I4_D[0] = 0.0;
    ppc->v_prob_I4_D[1] = 0.0;
    ppc->v_prob_I4_D[2] = 0.0001;
    ppc->v_prob_I4_D[3] = 0.0001;
    ppc->v_prob_I4_D[4] = 0.0005;
    ppc->v_prob_I4_D[5] = 0.003;
    ppc->v_prob_I4_D[6] = 0.01;
    ppc->v_prob_I4_D[7] = 0.02;
    ppc->v_prob_I4_D[8] = 0.04;

    // set the probability of progression from the "HA_i" -state to the CA state, for each HA_i stage individually
    double ph2c = 1.0 - pow( 1.0 - G_CLO_P_HOSP_TO_ICU, ( 1.0/((double)NUMHA) ) );
    ppc->v_prob_HA_CA[0] = 0.0;
    ppc->v_prob_HA_CA[1] = 0.0;     // no one under 20 progresses to the ICU from the Hospitalized State
    ppc->v_prob_HA_CA[2] = ph2c/10.0; 
    ppc->v_prob_HA_CA[3] = ph2c/10.0;
    ppc->v_prob_HA_CA[4] = ph2c/5.0;
    ppc->v_prob_HA_CA[5] = ph2c/2.0;
    ppc->v_prob_HA_CA[6] = ph2c;
    ppc->v_prob_HA_CA[7] = ph2c;
    ppc->v_prob_HA_CA[8] = ph2c;

    // set the probability of death for HA4 for all 9 age classes - 
    ppc->v_prob_HA4_D[0] = 0.0;
    ppc->v_prob_HA4_D[1] = 0.0;
    ppc->v_prob_HA4_D[2] = 0.0001;
    ppc->v_prob_HA4_D[3] = 0.0001;
    ppc->v_prob_HA4_D[4] = 0.0005;
    ppc->v_prob_HA4_D[5] = 0.004;
    ppc->v_prob_HA4_D[6] = 0.008;
    ppc->v_prob_HA4_D[7] = 0.03;
    ppc->v_prob_HA4_D[8] = 0.05;

    // set the probability of death for CA-individuals for all 9 age classes - 
    ppc->v_prob_CA_D[0] = 0.0; // these can be set to about .125 for the higher age classes; from the Seattle ICU paper on 24 patients
    ppc->v_prob_CA_D[1] = 0.001;
    ppc->v_prob_CA_D[2] = 0.125;
    ppc->v_prob_CA_D[3] = 0.125;
    ppc->v_prob_CA_D[4] = 0.125;
    ppc->v_prob_CA_D[5] = 0.125;
    ppc->v_prob_CA_D[6] = 0.125;
    ppc->v_prob_CA_D[7] = 0.125;
    ppc->v_prob_CA_D[8] = 0.125;
    // set the probability of death for CA-individuals for all 9 age classes - 
    ppc->v_prob_CA_V[0] = 0.0; // these can be set to about .75 for the higher age classes; from the Seattle ICU paper on 24 patients
    ppc->v_prob_CA_V[1] = 0.75;
    ppc->v_prob_CA_V[2] = 0.75;
    ppc->v_prob_CA_V[3] = 0.75;
    ppc->v_prob_CA_V[4] = 0.75;
    ppc->v_prob_CA_V[5] = 0.75;
    ppc->v_prob_CA_V[6] = 0.75;
    ppc->v_prob_CA_V[7] = 0.75;
    ppc->v_prob_CA_V[8] = 0.75;

    // from the Seattle ICU data on 24 patients, this probability is 60% -- obviously, it's a small sample size of older patients
    ppc->v_prob_V_D[0] = 0.0; 
    ppc->v_prob_V_D[1] = 0.001;
    ppc->v_prob_V_D[2] = 0.60;
    ppc->v_prob_V_D[3] = 0.60;
    ppc->v_prob_V_D[4] = 0.60;
    ppc->v_prob_V_D[5] = 0.60;
    ppc->v_prob_V_D[6] = 0.60;
    ppc->v_prob_V_D[7] = 0.60;
    ppc->v_prob_V_D[8] = 0.60;// TODO INCLUDE IN THE DIFF EQUATIONS

    
    for(int ac=0; ac<NUMAC; ac++) ppc->v_fraction_asymp[ac]=0.25;
    
    //printf("\n\n\t %1.3f \t %1.3f \t %1.3f \n\n", ppc->v[0],ppc->v[1],ppc->v[2]); fflush(stdout);

    
    
    //
    // ###  3.  RUN THE MODEL
    //
    generate_trajectories( 0.01, 3.0, 0.0, G_CLO_TF, 0.0 );


    if( OutFile != NULL ) fclose(OutFile);
    
    
    if( G_B_DIAGNOSTIC_MODE )
    {
        int ac;
        
        printf("\n\t\t\t\t 0-9 \t\t 10-19  \t 20-29  \t 30-39  \t 40-49  \t 50-59  \t 60-69  \t 70-79  \t 80+ ");
        printf("\n\tTotal Symp Cases");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTJ + ac]+0.5) );
        }
        printf("\n\tTotal Deaths    ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTD + ac]+0.5) );
        }
        printf("\n\tCFR             ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%1.3f%%   ", 100.0 * yic[STARTD + ac] / yic[STARTJ + ac]  );
        }

        printf("\n\n\tTotal Hosp Cases");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTK + ac]+0.5) );
        }
        printf("\n\tHosp FR         ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%1.1f%%     ", 100.0 * yic[STARTD + ac] / yic[STARTK + ac]  );
        }
        
    
        printf("\n\n");
    }
    

    delete[] yic;
    delete ppc;
    return 0;
}






