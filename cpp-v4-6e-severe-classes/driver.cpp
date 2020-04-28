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
double G_CLO_INTRODUCTION_TIME = -1;
int G_CLO_INTRODUCTION_COUNT;
string G_CLO_LOCATION = "RI";

double G_CLO_TF = 365.0;
double G_CLO_P_HOSP_TO_ICU = 0.30;

double G_CLO_SYMP_FRAC = 0.25;
double G_CLO_HOSPFRAC_YOUNG_DEV = 1.0;
double G_CLO_HOSPFRAC_MID_DEV = 0.5;
double G_CLO_HOSPFRAC_OLD_DEV = 0.5;

double G_CLO_ICUFRAC_DEV = 1.0;

double G_CLO_VENTDEATH_MID_DEV = 0.7;

double G_CLO_RELATIVE_BETA_HOSP = 0.2;

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
    ppc->v[ i_len_incub_period ]                            = 6.0;  // this is an average of the Lauer et al estimate (5.5d) and Backer et al (6.5d)
    ppc->v[ i_len_symptomatic_infectious_period_phase_1 ]   = 7.0;
    ppc->v[ i_len_symptomatic_infectious_period_phase_2 ]   = 7.0;
    ppc->v[ i_len_medicalfloor_hospital_stay ]              = 10.7; //   take from Lewnard et al, survivors only
    
    // params below are for relative infectiousness of certain individuals; I1 to I4 individuals have infectiousness = 1.0
    ppc->v[ i_phi_asymp ] = 0.5;        // set to same value as Imperial college models
    ppc->v[ i_phi_incub ] = 0.5;        // currently a complete unknown
    ppc->v[ i_phi_hosp ]  = 1.0;
    ppc->v[ i_phi_icu ]   = 1.0;
    ppc->v[ i_phi_vent ]  = 0.1;         // to discuss
    
    // params below are for relative contact levels of hospitalized, ICUed, and vented individuals
    ppc->v[ i_beta_hosp ] = G_CLO_RELATIVE_BETA_HOSP * ppc->v_betas[0];  // these relbeta's should do not change with social distancing measures
    ppc->v[ i_beta_icu ]  = G_CLO_RELATIVE_BETA_HOSP * ppc->v_betas[0];  // because the contact rate of a hospitalized patient is unaffected by the 
    ppc->v[ i_beta_vent ] = G_CLO_RELATIVE_BETA_HOSP * ppc->v_betas[0];  // social distancing policies set for health individuals
    
    
    //for(int ac=0; ac<NUMAC; ac++) ppc->v_fraction_asymp[ac]=0.25;

    // set the fraction of individuals who progress to asymptomatic infection
    double a=G_CLO_SYMP_FRAC;  //NOTE this number has to be somewhere between 0.1 and 0.3, closer to 0.3 probably; default is 0.25;
    //
    ppc->v_fraction_asymp[0] = 1.0 - a * 0.05; // 98.5% to 99.5% asymp       // these relative params are all taken from Joseph Wu et al, Nat Med, 2020
    ppc->v_fraction_asymp[1] = 1.0 - a * 0.08; // 97.6% to 99.2% asymp
    ppc->v_fraction_asymp[2] = 1.0 - a * 0.41; // 87.7% to 95.9% asymp
    ppc->v_fraction_asymp[3] = 1.0 - a * 1.00; // 70% to 90% asymp           // 30-39, the reference group
    ppc->v_fraction_asymp[4] = 1.0 - a * 1.35; // 59.5% to 86.5% asymp
    ppc->v_fraction_asymp[5] = 1.0 - a * 1.99; // 40.3% to 80.1% asymp 
    ppc->v_fraction_asymp[6] = 1.0 - a * 2.85; // 14.5% to 71.5% asymp 
    ppc->v_fraction_asymp[7] = 1.0 - a * 3.05; //  8.5% to 69.5% asymp 
    ppc->v_fraction_asymp[8] = 1.0 - a * 2.49; // 25.3% to 75.1% asymp 

    // set the fraction of individuals who are hospitalized immediately after I_2
    //
    double b1 = G_CLO_HOSPFRAC_YOUNG_DEV;   // default is 1.8 :: REASON is that we want these rates to match the hosp-age-dist in the Lewnard paper
    double b2 = G_CLO_HOSPFRAC_MID_DEV;     // default is 
    double b3 = G_CLO_HOSPFRAC_OLD_DEV;     // default is 
    ppc->v_fraction_hosp[0] = 0.025*b1;     // NOTE all of these numbers are taken directly from the CDC MMRW Mar 27 report
    ppc->v_fraction_hosp[1] = 0.025*b1;
    ppc->v_fraction_hosp[2] = 0.208*b2;
    ppc->v_fraction_hosp[3] = 0.208*b2;
    ppc->v_fraction_hosp[4] = 0.246*b2;
    ppc->v_fraction_hosp[5] = 0.292*b2;
    ppc->v_fraction_hosp[6] = 0.368*b3;
    ppc->v_fraction_hosp[7] = 0.511*b3;
    ppc->v_fraction_hosp[8] = 0.645*b3;

    // TODO deprecate this parameter - this will not be used
    // set the fraction of individuals who are immediately admitted to critical care after I_2
    ppc->v_fraction_crit[0] = 0.0;
    ppc->v_fraction_crit[1] = 0.0;
    ppc->v_fraction_crit[2] = 0.0;
    ppc->v_fraction_crit[3] = 0.0;
    ppc->v_fraction_crit[4] = 0.0;
    ppc->v_fraction_crit[5] = 0.0;
    ppc->v_fraction_crit[6] = 0.0;
    ppc->v_fraction_crit[7] = 0.0;
    ppc->v_fraction_crit[8] = 0.0;

    // set the probability of death for I4 for all 9 age classes
    ppc->v_prob_I4_D[0] = 0.0; // COMPLETELY UNKNOWN
    ppc->v_prob_I4_D[1] = 0.0;
    ppc->v_prob_I4_D[2] = 0.0;
    ppc->v_prob_I4_D[3] = 0.0;
    ppc->v_prob_I4_D[4] = 0.0;
    ppc->v_prob_I4_D[5] = 0.0;
    ppc->v_prob_I4_D[6] = 0.0;
    ppc->v_prob_I4_D[7] = 0.0;
    ppc->v_prob_I4_D[8] = 0.0;

    // set the probability of progression from the "HA_i" -state to the CA state, for each HA_i stage individually
    // NOTE THIS EQ USED TO CALCULATE NUMBERS BELOW: double ph2c = 1.0 - pow( 1.0 - G_CLO_P_HOSP_TO_ICU, ( 1.0/((double)NUMHA) ) );
    //          CALCULATIONS BELOW ASSUME THAT NUMHA=4
    ppc->v_prob_HA_CA[0] = 1.0 - pow( 1.0 - 0.304,  ( 1.0/((double)NUMHA) ) ); 
    ppc->v_prob_HA_CA[1] = 1.0 - pow( 1.0 - 0.293,  ( 1.0/((double)NUMHA) ) );    // NOTE the estimated Lewnard et al (medRxiv, Apr 16) probabilities are used here
    ppc->v_prob_HA_CA[2] = 1.0 - pow( 1.0 - 0.2825, ( 1.0/((double)NUMHA) ) );    //      averaged over male/female equally
    ppc->v_prob_HA_CA[3] = 1.0 - pow( 1.0 - 0.301,  ( 1.0/((double)NUMHA) ) );
    ppc->v_prob_HA_CA[4] = 1.0 - pow( 1.0 - 0.463,  ( 1.0/((double)NUMHA) ) );
    ppc->v_prob_HA_CA[5] = 1.0 - pow( 1.0 - 0.4245, ( 1.0/((double)NUMHA) ) );
    ppc->v_prob_HA_CA[6] = 1.0 - pow( 1.0 - 0.460,  ( 1.0/((double)NUMHA) ) );
    ppc->v_prob_HA_CA[7] = 1.0 - pow( 1.0 - 0.4835, ( 1.0/((double)NUMHA) ) );
    ppc->v_prob_HA_CA[8] = 1.0 - pow( 1.0 - 0.416,  ( 1.0/((double)NUMHA) ) );
    
    // the probabilities above range from: 0.07 to 0.16
    double c1 = G_CLO_ICUFRAC_DEV;
    for(int ac=0; ac<NUMAC; ac++)
    {
        ppc->v_prob_HA_CA[ac] *= c1;
        //ppc->v_prob_HA_CA[ac] = 0.0;
    }
    
    // set the probability of death for HA4 for all 9 age classes - 
    ppc->v_prob_HA4_D[0] = 0.0;
    ppc->v_prob_HA4_D[1] = 0.0;
    ppc->v_prob_HA4_D[2] = 0.0;
    ppc->v_prob_HA4_D[3] = 0.0;
    ppc->v_prob_HA4_D[4] = 0.0;
    ppc->v_prob_HA4_D[5] = 0.0;   // NOTE there are no data right now for these numbers
    ppc->v_prob_HA4_D[6] = 0.0;
    ppc->v_prob_HA4_D[7] = 0.03;
    ppc->v_prob_HA4_D[8] = 0.05;

    
    // set the probability of ventilation for CA-individuals for all 9 age classes - 
    ppc->v_prob_CA_V[0] = 0.75; // these can be set to about .75 for the higher age classes; from the Seattle ICU paper on 24 patients
    ppc->v_prob_CA_V[1] = 0.75;
    ppc->v_prob_CA_V[2] = 0.75;
    ppc->v_prob_CA_V[3] = 0.75;
    ppc->v_prob_CA_V[4] = 0.75;
    ppc->v_prob_CA_V[5] = 0.75;
    ppc->v_prob_CA_V[6] = 0.75;
    ppc->v_prob_CA_V[7] = 0.75;
    ppc->v_prob_CA_V[8] = 0.75;

    // from the Seattle ICU data on 24 patients, this probability is 60% -- obviously, it's a small sample size of older patients
    // these are being set to the ICU-to-Death probabilities since it's very difficult to get good data on death when on and not on a ventilator (for ICU patients)
    //
    double vd = G_CLO_VENTDEATH_MID_DEV; // default set to 1.0
    ppc->v_prob_V_D[0] = 0.03125;       // NOTE set directly from the Lewnard paper; very little data here
    ppc->v_prob_V_D[1] = 0.05119;       // NOTE set directly from the Lewnard paper; very little data here
    ppc->v_prob_V_D[2] = 0.15;          // this range should be between 14% (Lewnard) and 16.7%  (Graselli)
    ppc->v_prob_V_D[3] = 0.15;          // this range should be between 13% (Lewnard) and 17%  (Graselli), but Yang LRM observed 0%
    ppc->v_prob_V_D[4] = 0.400*vd;          // this range should be between 31% and 50%  (Lewnard, Graselli, Yang LRM)
    ppc->v_prob_V_D[5] = 0.460*vd;          // this range should be between 26% and 70%  (Lewnard, Graselli, Yang LRM)
    ppc->v_prob_V_D[6] = 0.585*vd;         // this range should be between 39% and 72%  (Lewnard, Graselli, Yang LRM, Bhatraju) 
    ppc->v_prob_V_D[7] = 0.70;          // this range should be between 60% and 88%  (Lewnard, Graselli, Yang LRM, Bhatraju)
    ppc->v_prob_V_D[8] = 0.90;          // this range should be between 60% and 100% (Lewnard, Graselli, Yang LRM, Bhatraju)


    // set the probability of death for CA-individuals for all 9 age classes 
    //
    // if you are in the ICU, and you don't progress to death, and your don't progress to ventilation, that means you are back on the medical-floor level of care
    //
    // very little data here; Seattle study on 24 ICU patients says this should be about 0.125
    //
    ppc->v_prob_CA_D[0] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];         
    ppc->v_prob_CA_D[1] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];                                     
    ppc->v_prob_CA_D[2] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];
    ppc->v_prob_CA_D[3] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];
    ppc->v_prob_CA_D[4] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];
    ppc->v_prob_CA_D[5] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];
    ppc->v_prob_CA_D[6] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];  // 0.146 --- so it's close to the 0.125 obersved in Seattle 24-patent study
    ppc->v_prob_CA_D[7] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];  // 0.175 --- so it's close to the 0.125 obersved in Seattle 24-patent study
    ppc->v_prob_CA_D[8] = (1.0-ppc->v_prob_CA_V[0]) * ppc->v_prob_V_D[0];  // 0.225 --- so it's not close to the average of 0.125 obersved in Seattle 24-patent study; but these patients are >80
    
    
    
    
    
    
    //
    // ###  3.  RUN THE MODEL
    //
    generate_trajectories( 0.01, 3.0, 0.0, G_CLO_TF, 0.0 );
    if( OutFile != NULL ) fclose(OutFile);




    
    
    //
    // ###  4.  OUTPUT DIAGNOSTICS
    //
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
            printf("\t %1.3f%%   ", 100.0 * yic[STARTD + ac] / yic[STARTJ + ac]  );
        }

        printf("\n\n\tTotal Hosp Cases");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t%7d ", (int)(yic[STARTK + ac]+0.5) );
        }
        printf("\n\tHosp FR         ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t  %1.1f%%   ", 100.0 * yic[STARTD + ac] / yic[STARTK + ac]  );
        }

        double total_num_hospitalized=0.0;
        for(ac=0; ac<NUMAC; ac++) total_num_hospitalized += yic[STARTK + ac];
        
        printf("\n\n\t%% Hosp by Age         ");
        for(ac=0; ac<NUMAC; ac++)
        {
            printf("\t  %1.1f%%   ", 100.0 * yic[STARTK + ac] / total_num_hospitalized  );
        }
        printf("\n\tLewnard et al  \t\t  0.1%% \t\t  0.2%% \t\t  3.6%% \t\t  9.8%% \t\t  14.6%% \t  20.5%% \t  22.2%% \t  16.6%% \t  12.3%%");
        
    
        printf("\n\n");
    }
    

    delete[] yic;
    delete ppc;
    return 0;
}






