//#include <iostream>
//#include <string>
//#include <cstdlib>

#include "essentials.h"
#include "assert.h"
#include "prms.h"


// constructor
prms::prms()
{
    v.insert( v.begin(), num_params, 0.0 );
    assert( v.size()==num_params );
    
    v_fraction_asymp.insert( v_fraction_asymp.begin(), NUMAC, 0.0 );
    v_fraction_hosp.insert(  v_fraction_hosp.begin(),  NUMAC, 0.0 );
    v_fraction_crit.insert(  v_fraction_crit.begin(),  NUMAC, 0.0 );
    v_prob_I4_D.insert(  v_prob_I4_D.begin(),  NUMAC, 0.0 );
    v_prob_HA4_D.insert( v_prob_HA4_D.begin(), NUMAC, 0.0 );
    v_prob_V_D.insert( v_prob_V_D.begin(), NUMAC, 0.0 );
    v_prob_CA_D.insert( v_prob_CA_D.begin(), NUMAC, 0.0 );
    v_prob_CA_V.insert( v_prob_CA_V.begin(), NUMAC, 0.0 );
        
    assert( v_fraction_asymp.size()==NUMAC );
    assert( v_fraction_hosp.size() ==NUMAC );
    assert( v_fraction_crit.size() ==NUMAC );
    assert( v_prob_I4_D.size() ==NUMAC );
    assert( v_prob_HA4_D.size()==NUMAC );
    assert( v_prob_V_D.size()==NUMAC );
    assert( v_prob_CA_D.size()==NUMAC );
    assert( v_prob_CA_V.size()==NUMAC );
    
    v_betas.clear();
    v_betatimes.clear();
    assert( v_betas.size() == 0 );
    assert( v_betatimes.size() == 0 );
    
    index_current_beta=-1;
    
}

// destructor
prms::~prms()
{
}


void prms::assign_new_beta( void )
{
    assert( v_betas.size() == v_betatimes.size() );
    assert( v_betas.size() > 0 );
    
    // if it's invalid, just set it to zero
    // this is what happens at initialization as well
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        index_current_beta = 0;
    }
    else // if it's valid, increment it if you can
    {
        if( index_current_beta != v_betas.size() - 1 ) // if it's NOT the last element, then increment it
        {
            index_current_beta++;
        }
    }
    
    // assign the possibly new beta value to the main parameter vector "v"
    v[i_beta] = v_betas[ index_current_beta ];
}


double prms::get_new_update_time( void )
{
    if( index_current_beta < 0 || index_current_beta >= v_betas.size() )
    {
        return 1000000.0;
    }
    else // if it's valid, increment it if you can
    {
        if( index_current_beta != v_betas.size() - 1 ) // if it's NOT the last element, then increment it
        {
            return v_betatimes[ index_current_beta+1 ];
        }
        else // if it is the last element in the array, just return 1000000.0;
        {
            return 1000000.0;
        }
    }
    
    
}


