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


##### -dev-hosp-young

This is the deviation in the parameter that determines the fraction of <60 individuals that are hospitalized. Default is 1.8.  Do not set this higher than 3.0. 

```diff
! Recommended prior distribution [0.5 - 2.5]
```

##### -dev-hosp-old

This is the deviation in the parameter that determines the fraction of >60 individuals that are hospitalized. Default is 1.0.  Do not set this higher than 1.3. 

```diff
! Recommended prior distribution [0.8 - 1.2]
```


##### -diag

This is a flag that prints out some basic diagnostic info as to who was infected, who was hospitalized, and who died.  To be used like this:

    `./odesim none -diag`

##### -symp-frac

This is the probability that a 30-39 year-old infected individual will develop symptoms. Default value is 0.25. If you set this higher than 0.325, it will be reset back down to 0.325.

```diff
! Recommended prior distribution [0.1 - 0.3]
```

##### -tf
Time at which the ODEs are stopped. Day 1 is Jan 1 2020.  So, if you want to run the simulation through to April 30 2020, you would call

   `./odesim none -tf 121`
   
The 'time initial' right now is set to zero, which means that it is set to Jan 1 2020.  Cases are introduced at a fixed time to begin the epidemic.  For now, we can use this fixed time, and I will add a command-line option later to change it.
