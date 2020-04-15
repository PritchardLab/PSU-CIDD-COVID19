### File describing the command-line usage of the executable 'odesim'

##### First argument is always the output filename

You can call `./odesim myoutput.tdl`

or `./odesim none`

and the 'none' keyword will skip outputting to any file.


##### -beta

You can call the function with an arbitrary number of beta parameters (transmission parameters) that allow for time-varying transmission in the model.  For example, you can call

   `./odesim none -beta 1.0 0.9 0.8 0.7`
   
and this will run the model with four different transmission rates, evenly distributed, in the order above, for the period between March 1 and the end of the simulation which is denoted by `tf = time final`.  


> WARNING
>
> For now, you have to set time final on the command line before you assign the betas
>
> e.g. `./odesim none -tf 120 -beta 3 4 5 6`

##### -diag

This is a flag that prints out some basic diagnostic info as to who was infecte, who was hospitalized, and who died.  To be used like this:

    `./odesim none -diag`

##### -ph2c

This is the probability that a hospitalized patient (in the acute stage of hospitalization) will at some point in the future need critical care (i.e. ICU care).


##### -tf
Time at which the ODEs are stopped. Day 1 is Jan 1 2020.  So, if you want to run the simulation through to April 30 2020, you would call

   `./odesim none -tf 121`
   
The 'time initial' right now is set to zero, which means that it is set to Jan 1 2020.  Cases are introduced at a fixed time to begin the epidemic.  For now, we can use this fixed time, and I will add a command-line option later to change it.
