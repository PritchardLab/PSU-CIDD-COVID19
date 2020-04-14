#include <math.h>
#include <stdio.h>
#include <ctype.h>
#include <vector>
#include <cstring>
#include <sstream>
#include <fstream>
#include <string>
//#include <assert.h>
#include <stdlib.h>
#include "essentials.h"
#include "prms.h"
#include "assert.h"



// BEGIN ### ### GLOBAL VARIABLES ### ###

extern double* yic;  
extern prms* ppc;
extern double G_CLO_INTRODUCTION_TIME;
extern int G_CLO_INTRODUCTION_COUNT;
extern FILE* OutFile;
extern double G_CLO_TF;
extern string G_CLO_LOCATION;

//  END  ### ### GLOBAL VARIABLES ### ###



bool isFloat( string myString );
void PrintUsageModes();
void SetLocationData( string loc );

bool isFloat( string myString ) {
    std::istringstream iss(myString);
    float f;
    iss >> noskipws >> f; // noskipws considers leading whitespace invalid
    // Check the entire string was consumed and if either failbit or badbit is set
    return iss.eof() && !iss.fail(); 
}

void PrintUsageModes()
{
    printf("\n\tUSAGE: ./odesim  outfilename   [-options] \n\n");
}


//
// parses command line arguments
void ParseArgs(int argc, char **argv)
{
    int i, start;
    start=2;

    if( argc<2 )
    { 
        PrintUsageModes(); 
	exit(-1);
    }
	
    if( argv[1][0]=='n' && argv[1][1]=='o' && argv[1][2]=='n' && argv[1][3]=='e' && argv[1][4]==0 )
    {
        //fprintf(stderr, "\n\n\tnot printing to outfile\n\n");
    }
    else 
    {
	OutFile = fopen( argv[1], "w" );
    }

    string str;

    i=start;
    
    // read in options
    while(i<argc)
    {
	str = argv[i];
        
        // ### 1 ### IF BLOCK FOR BETA
        if( str == "-beta" )
        {
            ppc->v_betas.clear();
            ppc->v_betatimes.clear();
            i++;
            bool bFirstBetaPushedBack = false;
            
            //BEGIN LOOPING THROUGH THE LIST OF BETAS
            while(i<argc)
            {
                string s( argv[i] );
                if( isFloat(s) ) // if the current string is a floating point number, write it into the v_betas array
                {
                    // if the command line argument is <0, just set it back to zero
                    double d = atof( argv[i] );
                    if( d < 0.0 ) d = 0.0;
                    
                    ppc->v_betas.push_back( d );
                    
                    // if it's the first beta value being read in, then push it back twice, as it will be used as 
                    // the Jan 1 to March 1 beta (when there were almost no cases) as well as the beta for the 
                    // first time period we are invesigating
                    if(!bFirstBetaPushedBack) {ppc->v_betas.push_back( d );bFirstBetaPushedBack=true;}
                    
                    // increment and move on in this sub-loop
                    i++;
                }
                else
                {
                    // if the current string is NOT a float, set the counter back down (so the string can be handled by the next if block
                    // and break out of the loop
                    i--;
                    break;
                }
                 
            } 
            //END OF LOOPING THROUGH THE LIST OF BETAS
             
            // MAKE SURE AT LEAST ONE BETA VALUE WAS READ IN 
            if( ppc->v_betas.size() == 0 )
            {
                fprintf(stderr,"\n\n\tIn parsing command-line arguments, '-beta' option must be followed by at least one floating point number.\n\n");
            }
            else
            {
                double total_time_period = G_CLO_TF - 60.0;
                double time_step = total_time_period / ( (double) ppc->v_betas.size() - 1 );
                
                ppc->v_betatimes.push_back(0.0);
                
                for(int jj=1; jj<ppc->v_betas.size(); jj++)
                {
                    ppc->v_betatimes.push_back( 60.0 + time_step*((double)(jj-1)) );
                }
                
            }
                         
        }
        // ### 2 ### IF BLOCK FOR TIME-FINAL    
        else if( str == "-tf" )
        {
            G_CLO_TF = atof( argv[++i] );
        }
        // ### 3 ### IF BLOCK FOR LOCATION
        else if( str == "-loc" )
        {
            G_CLO_LOCATION = argv[++i];
        }
        // ### FINAL ### IF BLOCK FOR AN UNKNOWN COMMAND-LINE OPTIONS            
 	else
 	{
 	    fprintf(stderr, "\n\tUnknown option [%s] on command line.\n\n", argv[i]);
 	    exit(-1);
 	}
        
            
        
        
        
        
        //END OF MAIN WHILE-LOOP BLOCK; INCREMENT AND MOVE ON
        i++;
    }
    
    // after parsing all of the arguments, do the following
    
    // 1. set the specific age structure, population, and introduction time for this location
    SetLocationData( G_CLO_LOCATION );


    return;
}

void SetLocationData( string loc )
{
    for(int i=0;i<STARTJ+NUMAC;i++) yic[i]=0.0; // zero everything out
    
    if( loc=="RI" )
    {
        // 2019 est pop is 1,059,361
        

        // proportions from https://www.statista.com/statistics/1022746/rhode-island-population-share-age-group/
        ppc->v[i_N] = 1059361.0;
        yic[0]  = ppc->v[i_N]   *   .105;   //  0-9
        yic[1]  = ppc->v[i_N]   *   .123;   //  10-19
        yic[2]  = ppc->v[i_N]   *   .140;   //  20-29
        yic[3]  = ppc->v[i_N]   *   .127;   //  30-39
        yic[4]  = ppc->v[i_N]   *   .124;   //  40-49
        yic[5]  = ppc->v[i_N]   *   .135;   //  50-59
        yic[6]  = ppc->v[i_N]   *   .120;   //  60-69
        yic[7]  = ppc->v[i_N]   *   .074;   //  70-79
        yic[8]  = ppc->v[i_N]-yic[0]-yic[1]-yic[2]-yic[3]-yic[4]-yic[5]-yic[6]-yic[7];   //  80+
        assert( yic[8] > 0.0 );
            
        G_CLO_INTRODUCTION_TIME = 65.0;
        G_CLO_INTRODUCTION_COUNT = 2;
        
    }
    else if( loc=="PA" )
    {
        
        
    }
    else
    {
        fprintf(stderr, "\n\tUnknown location [%s] entered after -loc on command line.\n\n", loc.c_str() );
        exit(-1);
    }
    
}





