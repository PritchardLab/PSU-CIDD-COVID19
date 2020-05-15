#ifndef PRMS
#define PRMS

#include <vector>
#include <math.h>


using namespace std;

// the i_phi parameters are the relative infectiousness of individuals in that state; the infectiousness of the I1 to I4 classes is assumed to be 1.0
//

enum parameter_index {i_N, i_beta, i_startday_beta, i_endday_beta, i_beta_hosp, i_beta_icu, i_beta_vent, i_len_incub_period, i_len_symptomatic_infectious_period_phase_1, i_len_symptomatic_infectious_period_phase_2, i_len_medicalfloor_hospital_stay, i_phi_asymp, i_phi_incub, i_phi_hosp, i_phi_icu, i_phi_vent, num_params}; 

typedef enum parameter_index prm_index;


class prms
{   
public:    
    explicit prms();          	// constructor
    ~prms();         	// destructor


    vector<double> v;			// this holds all the parameters; they are indexed by the enums above

    //BEGIN --- FOR ALL THE VECTORS BELOW, THERE WILL BE NUMAC ENTRIES FOR THE NINE AGE CLASSES
    //
    
    vector<double> v_fraction_asymp;    
    
    vector<double> v_fraction_hosp;     // fraction of infected & symp individuals who are hospitalized after class I_2
    vector<double> v_fraction_crit;     // fraction of infected & symp individuals who are immediately admitted to critical care after class I_2
    
    vector<double> v_prob_I4_D;         // probability of death for a non-hospitalized person
    
    vector<double> v_prob_HA4_D;        // probability of death on the regular medical-level-of-care hospital floor
    vector<double> v_prob_HA_CA;        // probability of progressing to the ICU from any of the HA-classes during the acute-phase hospitalization phase
    
    vector<double> v_prob_V_D;          // probability of death if you are on a ventilator

    
    vector<double> v_prob_CA_D;         // probability of death directly from the CA-class (critical care, acute phase, but not on ventilator)
    vector<double> v_prob_CA_V;         // probability of progressing to mechanical ventilation from the CA-class 
                                        // if you take 'one minus' the two probabilities above, you will get the probability that someone moves from CA to HR
    
    

    //
    //END --- FOR ALL THE VECTORS ABOVE, THERE ARE NUMAC ENTRIES FOR THE NINE AGE CLASSES
    
    
    // NOTE THE TWO ARRAYS BELOW MUST BE OF THE SAME LENGTH
    //      THE "BETATIME" IS THE START-TIME THAT THE CORRESPONDING BETA COMES INTO EFFECT
    //      THE LAST BETA STAYS IN EFFECT FOREVER ONCE WE HAVE PASSED THE LAST "BETA-TIME"
    vector<double> v_betas;             // these are all the different beta parameters that will change every T time steps
                                        // where the goal is to have T=1 so you have beta changing daily
                                        
    vector<double> v_betatimes;         // these are all the different times that the beta parameters change
    
    
    
    
    int index_current_beta;
    void assign_new_beta( void );
    double get_new_update_time( void );
    
    
};

#endif // PRMS
