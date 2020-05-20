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

This is the deviation in the parameter that determines the fraction of <20 individuals that are hospitalized. Default is 1.8.  Do not set this higher than 3.0. 

```diff
! Recommended prior distribution [0.5 - 2.5]
```

##### -dev-hosp-mid

This is the deviation in the parameter that determines the fraction of 20-60 individuals that are hospitalized. Default is 1.0.  Do not set this higher than 3.0. 

```diff
! Recommended prior distribution [0.1 - 2.0]
```

##### -dev-hosp-old

This is the deviation in the parameter that determines the fraction of >60 cases that are hospitalized. Despite the fact that there are a lot of data to parametrize the hospitalization rates for older individuals, this parameter is not easily identifiable with the reporting rate. This means that hospitalization fractions, in some contexts, may appear low due to outreach and testing in nursing homes or widespread random testing in the community. Default is 0.5.  Do not set this higher than 1.3. 

```diff
! Recommended prior distribution [0.1 - 1.2]
```


##### -dev-icu-frac

This is the deviation in the parameter that determines the fraction of hospitalized individuals that are admitted to the ICU. This deviation will apply across all age groups.  Default is 1.0.  If you set this higher than 5.0, it will be reset back down to 5.0 

```diff
! Recommended prior distribution [0.5 - 1.5]
```

##### -dev-len-hospstay

This is the deviation in the parameter that determines the length of the typical hospital stay on the medical floor.  Default is 1.0.  If you set this higher than 2.5, it will be reset back down to 2.5 

```diff
! Recommended prior distribution [0.5 - 2.0]
```


##### -dev-ventdeath-mid

This is the deviation in the parameter that determines the fraction of 40-70 individuals that die after being on a ventilator. Default is 0.7.  Do not set this higher than 1.7. 

```diff
! Recommended prior distribution [0.5 - 1.3]
```


##### -diag

This is a flag that prints out some basic diagnostic info as to who was infected, who was hospitalized, and who died.  To be used like this:

    `./odesim none -diag`

##### -introday

The introduction day of the first case. Feb 1 2020 is day 32, March 1 2020 is day 61, etc.


##### -loc

The location currently being modeled. E.g.

    `./odesim none -loc PA`
    
for Pennsylvania. Acceptable arguments currently are PA, MA, and RI.  Default is RI.


##### -rel-beta-hosp

This is the relative beta, i.e. the relative population mixing paramter, for a hospitalized individual (including ICU and ventilated). Default is 0.2. Must be between 0.0 and 1.0. Note that this beta does not change when social distancing is enforced.

```diff
! Recommended prior distribution [0.0 - 0.3]
```

##### -symp-frac

This is the probability that a 30-39 year-old infected individual will develop symptoms. Default value is 0.25. If you set this higher than 0.325, it will be reset back down to 0.325.

```diff
! Recommended prior distribution [0.1 - 0.3]
```

##### -tf
Time at which the ODEs are stopped. Day 1 is Jan 1 2020.  So, if you want to run the simulation through to April 30 2020, you would call

   `./odesim none -tf 121`
   
The 'time initial' right now is set to zero, which means that it is set to Jan 1 2020.  Cases are introduced via the `-introday` command-line argument.  
