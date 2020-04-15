#ifndef ESSENTIALS
#define ESSENTIALS


#include <math.h>
#include <stdio.h>

// This is the number of age classes in the model, 9 age classes currently: 0-9,10-19, etc up to 80+
#define NUMAC 9


// number of stages for the exposed class (E-class) - these are individuals who have recently become exposed (and thus infected)
// but they are still in their incubation period, they are not showing symptoms, and they are likely not infectious,
// except maybe in the late stages of the incubation period
//
// it is currently set to 6 to give the incubation period a coefficiant of variation of 1/sqrt(6) which is about 0.4
#define NUME 6
#define STARTE 9    // this is the start-index of the exposed (E) individuals;
                    // this has to be NUMAC, because you have the 9 S-classes in the index before you start indexing the E-classes

                    
// number of stages for the asymptomatic class (A-class); these are individuals who are either (1) fully asymptomatic, or (2) have symptoms that
// are mild enough that they would not report to a doctor/clinic under normal circumstances
#define NUMA 4
#define STARTA 63   // this is the start-index of for the asymptomatic individuals


// number of infected classes (I-classes) - these individuals are symptomatic and infectious (when you transition from E6 to I1, that
// transition coincides with the arrival of symptoms
// - these individuals are not hospitalized
// - this class includes individuals that have no dyspnea, no hypoxia, and no other underlying conditions that would cause a hospitalization
// - this class includes individuals that have dyspnea, no hypoxia, and no other underlying conditions that would cause a hospitalization; these individuals would be 
//   treated as out-patient at a hospital or clinic
#define NUMI 4
#define STARTI 99   // this is the start-index of for the infected individuals
                    // this has to be NUMAC + NUMAC*NUME + NUMAC*NUMA because the order of the indexing goes S, then E, then A, then I

                    
// number of Hospitalized individuals in Acute phase (i.e. not recovering and getting better, but just arrived and getting worse)
#define NUMHA 4
#define STARTHA 135


// number of individuals in Critical Care (ICU) in acute phase (i.e. not recovering)
#define STARTCA 171                    
                    

// V = individuals on mechanical ventilation; NUMV is the number of stages in class V
#define NUMV 6
#define STARTV 180


// number of individuals in Critical Care (ICU) in recovering phase (post ventilation)
#define STARTCR 234                    


// number of hospitalized individuals in recovering phase (post ICU)
#define STARTHR 243                    


// D = dead individuals; there will be NUMAC dead classes.  No need to split the dead individuals into separate stages.  
#define STARTD 252
                    
                    
// there will be NUMAC recovered classes.  No need to split the recovered into separate groups
#define STARTR 261


// ### CLASSES THAT KEEP TRACK OF CUMULATIVE SYMPTOMATIC INCIDENCE
// number of J-classes - these are the cumulative symptomatic incidence classes
#define STARTJ 270

// ### CLASSES THAT KEEP TRACK OF CUMULATIVE HOSPITALIZATION INCIDENCE
// these are the K-classes 
#define STARTK 279


// number of equations 9 S classes, 54 E-classes, 36 A-classes, 36 I-classes, 9 R-classes, and 9 J-classes
//
// NOTE NOTE NOTE TOTAL NUMBER OF EQUATIONS IS 288






//#define HAVE_INLINE 	// this will make GSL vector operations faster

using namespace std;

//
//
// this function contains the ode system
//
int func(double t, const double y[], double f[], void *params);



//void* jac;	// do this for C-compilation
//
// for C++ compilation we are replacing the void* declaration above with
// the inline dummy declaration below
inline int jac(double a1, const double* a2, double* a3, double* a4, void* a5)
{
    return 0;	
};



inline double agesum( double* p, int numac=NUMAC )
{
    double sum=0;
    for( int i=0; i < numac; i++ ) sum += p[i];
    return sum;
}






#endif // ESSENTIALS
