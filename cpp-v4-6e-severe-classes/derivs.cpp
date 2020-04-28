#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "derivs.h"
#include "essentials.h"
#include "prms.h"

extern prms* ppc;



void derivs( double t, double *y, double *dydt)
{
    int i,ac;

    //printf("\n\n\t ************* %1.3f \t %1.3f \t %1.3f \n\n", ppc->v[0],ppc->v[1],ppc->v[2]); fflush(stdout);

    
    // ### 0 ### compute the force of infection 
    //
    double foi=0.0;

    for(i=NUME+4; i<NUME; i++) // loop through the last two stages of exposed individuals
    {
        for(ac=0; ac<NUMAC; ac++)
        {
            foi += ppc->v[i_beta] * ppc->v[i_phi_incub] * y[STARTE + i*NUMAC + ac];   // STARTI is the starting index of all of the I-classes 
        }
    }
    for(i=0; i<NUMAC*NUMA; i++) // loop through all asymptomatic individuals
    {
        foi += ppc->v[i_beta] * ppc->v[i_phi_asymp] * y[STARTA + i];   // STARTI is the starting index of all of the I-classes 
    }
    for(i=0; i<NUMAC*NUMI; i++) // you want to loop across all NUMI stages and all NUMAC ages of infected individuals
    {
        foi += ppc->v[i_beta] * y[STARTI + i];      // STARTI is the starting index of all of the I-classes 
    }                                               // no need to multiply by a relative infectiousness parameter because this is 1.0 here (the reference case)
    for(i=0; i<NUMAC*NUMHA; i++) // loop through all hospitalized individuals
    {
        foi += ppc->v[i_beta_hosp] * ppc->v[i_phi_hosp] * y[STARTHA + i];   
    }
    for(i=0; i<NUMAC; i++) // loop through all ICU individuals
    {
        foi += ppc->v[i_beta_icu] * ppc->v[i_phi_icu] * y[STARTCA + i];   
    }
    for(i=0; i<NUMV*NUMAC; i++) // loop through all ventilated individuals
    {
        foi += ppc->v[i_beta_vent] * ppc->v[i_phi_vent] * y[STARTV + i];   
    }
    for(i=0; i<NUMAC; i++) // loop through all CR individuals
    {
        foi += ppc->v[i_beta_icu] * ppc->v[i_phi_icu] * y[STARTCR + i];   
    }
    for(i=0; i<NUMAC; i++) // loop through all HR individuals
    {
        // foi += ppc->v[i_beta_hosp] * ppc->v[i_phi_hosp] * y[STARTHR + i];    // NOTE this is set to zero-contribution to the FOI right now because        
    }                                                                           // these individuals are on day 10-15 of their infection, and Wolfel et al (Nature, 2020) suggest that
                                                                                // positivity by virus-culture should be low by this time 
    
    
    
    
    foi *= ppc->v[i_beta]; // this is the beta parameters
    
    double popsize=ppc->v[i_N];

    
    //NOTE the variables below are transitions rate (tr) variables between different compartment types; e.g. trv is the transition rate among V-classes
    //
    // this is the transition rate among the E-classes 
    double tre = ((double)NUME) / ppc->v[i_len_incub_period];
    
    // this is the transition rate among the I-classes
    double tri1 = (((double)NUMI)/2.0) / ppc->v[ i_len_symptomatic_infectious_period_phase_1 ];
    double tri2 = (((double)NUMI)/2.0) / ppc->v[ i_len_symptomatic_infectious_period_phase_2 ];
    
    // this is the transition rate among the HA-classes
    double trha = ((double)NUMHA) / ppc->v[ i_len_medicalfloor_hospital_stay ];
    
    
    // this is the transition rate among the A-classes  TODO update this; it is hardcoded as 6 days for now
    //double tri = 0.5;
    double tra = 0.75;
    
    // this is the transition rate out of the 'acute ICU' stage meaning your mean time here is two days
    double trca = 0.5;
    
    // this is the transition rate among the V-classes ... each stage is about 1.8 days
    // meaning you have about 10.8 days on a ventilator in the 6 V classes
    double trv = 1.0/1.8;
    
    double trhr = 0.4;  //NOTE BethG says this should be about 2-3 days.
    double trcr = 0.50; //TODO this is a placeholder; assign this somewhere
    
    
    
        
    // ### 1 ### from index=0 to index=NUMAC-1 you have the S-classes, one for each age class
    //
    // first the S-classes
    for(ac=0;ac<NUMAC;ac++)
    {
        dydt[ac] = - foi * y[ac] / popsize; 
    }

    // ### 2 ### from index=NUMAC to index=NUMAC+NUME*NUMAC you have the E-classes, one for each age class
    // there should be 54 E-classes (54 = NUME*NUMAC)
    for(i=0;i<NUME;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the E_1 class and the S-classes have to flow into it
            {
                dydt[STARTE + i*NUMAC + ac] = foi * y[ac] / popsize - tre * y[STARTE + i*NUMAC + ac];
            }
            else // these are the classes E_2 and higher
            {
                dydt[STARTE + i*NUMAC + ac] = tre * y[STARTE + (i-1)*NUMAC + ac] - tre * y[STARTE + i*NUMAC + ac];
            }
            
        }
    }

    
    // ### 3 ###    ASYMPTOMATIC, INFECTED, SOMEWHAT INFECTIOUS CLASSES (the A-class)
    //              there should be NUMAC*NUMA = 36 of these as well 
    for(i=0;i<NUMA;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the A_1 class and a fraction of the E_final class has to flow into it
            {
                dydt[STARTA + i*NUMAC + ac] = ppc->v_fraction_asymp[ac]*tre*y[STARTE + (NUME-1)*NUMAC + ac] - tra * y[STARTA + i*NUMAC + ac]; 
            }
            else // meaning this is the A_2 class or higher
            {
                dydt[STARTA + i*NUMAC + ac] = tra*y[STARTA + (i-1)*NUMAC + ac] - tra * y[STARTA + i*NUMAC + ac];
            }
            
        }
    }
    
    
    
    // ### 4 ###    INFECTED, INFECTIOUS, AND SYMPTOMATIC CLASSES (THE I-CLASS)
    // there are 36 I-classes (36 = NUMI*NUMAC)
    for(i=0;i<NUMI;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the I_1 class and a fraction of the E_final class has to flow into it
            {
                dydt[STARTI + i*NUMAC + ac] = (1.0-ppc->v_fraction_asymp[ac])*tre*y[STARTE + (NUME-1)*NUMAC + ac] - tri1 * y[STARTI + i*NUMAC + ac]; 
            }
            else if(i==1) // meaning this is the I_2 class 
            {
                dydt[STARTI + i*NUMAC + ac  ] = tri1 * y[STARTI + (i-1)*NUMAC + ac] - tri1 * y[STARTI + i*NUMAC + ac];
            }
            else if(i==2) // meaning this is the I_3 class .... NOTE you have to be **non-hospitalized** to flow into this class
            {
                dydt[STARTI + i*NUMAC + ac] = (1.0-ppc->v_fraction_hosp[ac]-ppc->v_fraction_crit[ac])*tri1 * y[STARTI + (i-1)*NUMAC + ac] - tri2 * y[STARTI + i*NUMAC + ac];
            }
            else if(i==3) // meaning this is the I_4 class 
            {
                dydt[STARTI + i*NUMAC + ac] = tri2 * y[STARTI + (i-1)*NUMAC + ac] - tri2 * y[STARTI + i*NUMAC + ac];
            }
            else
            {
                assert(false);   
            }
                        
        }
    }

    
    
    // ### 5 ###    HOSPITALIZED INDIVIDUALS, ACUTE-PHASE, MEANING THAT THEY RECENTLY ARRIVED AT THE HOSPITAL AND THEY ARE AT RISK FOR WORSENING CONDITIONS
    //              there should be NUMAC*NUMHA = 36 of these as well 
    for(i=0;i<NUMHA;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the HA_1 class and a fraction of the I_2 class flows into here
            {
                dydt[STARTHA + i*NUMAC + ac] = ppc->v_fraction_hosp[ac]*tri1*y[STARTI + NUMAC + ac] - trha * y[STARTHA + i*NUMAC + ac]; 
            }
            else // meaning this is the HA_2 class or higher
            {
                dydt[STARTHA + i*NUMAC + ac] = (1.0 - ppc->v_prob_HA_CA[ac]) * trha * y[STARTHA + (i-1)*NUMAC + ac] - trha * y[STARTHA + i*NUMAC + ac];
            }
            
        }
    }
    
    
    
    // ### 6 ###    CRITICAL-CARE/ICU INDIVIDUALS, ACUTE-PHASE, MEANING THAT THEY RECENTLY ARRIVED IN THE ICU AND THEY ARE AT RISK FOR WORSENING CONDITIONS
    //              THERE IS ONLY ONE "STAGE" OF ACUTE ICU FOR NOW
    for(i=0;i<1;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            dydt[STARTCA + ac] = 0.0;
            
            for(int stg=0; stg<NUMHA; stg++)
            {
                dydt[STARTCA + ac] += ppc->v_prob_HA_CA[ac] * trha * y[STARTHA + stg*NUMAC + ac];
            }
            
            dydt[STARTCA + ac] += ppc->v_fraction_crit[ac]*tri1*y[STARTI + NUMAC + ac]; //NOTE this should be zero; remember to deprecate this
            
            dydt[STARTCA + ac] -= trca * y[STARTCA + ac];            
        }
    }
    
    

    
    // ### 7 ###    INDIVIDUALS ON VENTILATORS
    //              there should be NUMAC*NUMV = 54 of these classes
    
    // this is the index of the 'special' V-class (currently, the fourth class) from which a patient can die
    int VK = 3;
    
    for(i=0;i<NUMV;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            if(i==0) // meaning this is the V_1 class and a fraction of the CA class flows into here
            {
                dydt[STARTV + i*NUMAC + ac] = ppc->v_prob_CA_V[ac]*trca*y[STARTCA + ac] - trv * y[STARTV + i*NUMAC + ac]; 
            }
            else if(i==VK+1) // meaning, this is the fifth class and patients who progress to death will not enter this class
            {
                dydt[STARTV + i*NUMAC + ac] = trv * (1.0-ppc->v_prob_V_D[ac]) * y[STARTV + (i-1)*NUMAC + ac] - trv * y[STARTV + i*NUMAC + ac];
            }
            else // meaning this is the V_2 class or higher
            {
                dydt[STARTV + i*NUMAC + ac] = trv * y[STARTV + (i-1)*NUMAC + ac] - trv * y[STARTV + i*NUMAC + ac];
            }
            
        }
    }

    
    // ### 8 ###    CRITICAL-CARE/ICU INDIVIDUALS, RECOVERING-PHASE, MEANING THAT THEY JUST GOT OFF A VENTILATOR
    //              THERE IS ONLY ONE "STAGE" OF RECOVERING ICU FOR NOW
    for(i=0;i<1;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            dydt[STARTCR + ac] = trv * y[STARTV + (NUMV-1)*NUMAC + ac] - trcr * y[STARTCR + ac];            
        }
    }

    
    // ### 9 ###    HOSPITALIZED INDIVIDUALS, RECOVERING-PHASE, MEANING THAT THEY JUST GOT HERE FROM THE ICU
    //              THERE IS ONLY ONE "STAGE" OF "HR" FOR NOW
    for(i=0;i<1;i++)
    {
        for(ac=0;ac<NUMAC;ac++)
        {
            dydt[STARTHR + ac] = trcr * y[STARTCR + ac]  + (1.0-ppc->v_prob_CA_V[ac]-ppc->v_prob_CA_D[ac]) * trca * y[STARTCA + ac] - trhr * y[STARTHR + ac]; 
            //TODO need to account for death here, some CRs won't enter here but will die
        }
    }
    
    
    
    // ### 10 ###    DEAD INDIVIDUALS (D) TODO need to allow for death of CR individuals
    //
    for(ac=0;ac<NUMAC;ac++)
    {
        // set deriv equal to zero
        dydt[STARTD + ac] = 0.0;
        
        // add asymptomatics recovering from A4 - WE ASSUME THAT THE DEATH RATE FORM THE ASYMPTOMATIC CLASS IS ZERO
        // dydt[STARTD + ac] += tra * y[STARTA + (NUMA-1)*NUMAC + ac]; 
        
        // add infecteds coming from I4 
        dydt[STARTD + ac] += ppc->v_prob_I4_D[ac] * tri2 * y[STARTI + (NUMI-1)*NUMAC + ac]; 

        // add HA-individuals coming from HA4 
        dydt[STARTD + ac] += ppc->v_prob_HA4_D[ac] * trha * y[STARTHA + (NUMHA-1)*NUMAC + ac]; 
        // add HR-individuals coming from HR 
        dydt[STARTD + ac] += ppc->v_prob_HA4_D[ac] * trhr * y[STARTHR + ac]; 

        // add CA-individuals coming from CA 
        dydt[STARTD + ac] += ppc->v_prob_CA_D[ac] * trca * y[STARTCA + ac]; 
        
        // add in V-individuals who will die (just V4 for now)
        dydt[STARTD + ac] += trv * ppc->v_prob_V_D[ac] * y[STARTV + VK*NUMAC + ac];
    }

    
    
    
    // ### 11 ###    RECOVERED INDIVIDUALS (R)
    //
    for(ac=0;ac<NUMAC;ac++)
    {
        // set deriv equal to zero
        dydt[STARTR + ac] = 0.0;
        
        // add asymptomatics recovering from A4
        dydt[STARTR + ac] += tra * y[STARTA + (NUMA-1)*NUMAC + ac]; 
        
        // add infecteds coming from I4
        dydt[STARTR + ac] += (1.0-ppc->v_prob_I4_D[ac]) * tri2 * y[STARTI + (NUMI-1)*NUMAC + ac]; 

        // add HA-individuals coming from HA4
        dydt[STARTR + ac] += (1.0 - ppc->v_prob_HA4_D[ac] - ppc->v_prob_HA_CA[ac]) * trha * y[STARTHA + (NUMHA-1)*NUMAC + ac];
        // add HR-individuals coming from HR
        dydt[STARTR + ac] += (1.0-ppc->v_prob_HA4_D[ac]) * trhr * y[STARTHR + ac];
        
    }
    
    
    
    // ### 12 ###    CUMULATIVE SYMPTOMATIC INCIDENCE CLASSES
    //               the 9 J-classes
    for(ac=0;ac<NUMAC;ac++)
    {
        dydt[STARTJ+ac] = (1.0-ppc->v_fraction_asymp[ac]) * tre * y[STARTE + (NUME-1)*NUMAC + ac];
    }
    
    
    // ### 13 ###    CUMULATIVE SYMPTOMATIC INCIDENCE CLASSES
    //               finally the 9 K-classes
    for(ac=0;ac<NUMAC;ac++)
    {
        dydt[STARTK+ac] = ( ppc->v_fraction_hosp[ac] + ppc->v_fraction_crit[ac] ) * tri1 * y[STARTI + NUMAC + ac] ;
    }
    
    
    
}

