[Work in progress]

Description
-----------

This repository is a CDCL version of [Mistral-2.0] dedicated to resolve (Disjunctive) Scheduling Problems. So far, only the standard Job-Shop problem is being (experimentally) tested. 

Installation
--------------
To use the Scheduler module, first compile the code using 

```sh
 make scheduler 
```


That's all! Try it out with 

```sh
bin/scheduler BENCHNAME -type jla
```

Where BENCHNAME is the instance file location and "-type jla" indicates that the file being used is in the format of Lawrence (the parameter should be omited with Taillard instances as the solver treats them by default). 
 
By default, a pure CP model is launched. If you want to enable Clause Learning, then a large number of parameters should be specified (although one is sufficient). So far, I'm using the following command line 

```sh
bin/scheduler BENCHNAME -type jla -fdlearning 2 -semantic 1 -keeplearning 1 -maxnogoodsize 12 -fixedLearntSize 1000  -fixedForget  5000  -fixedlimitSize   50000 -probforget 90 -vsids 1 -forgetdecsize  8
```


Parameters
----------

Several parameters can be used. These are probably the most important ones used with Clause Learning (I'll update the list later) : 

* -fdlearning : Used to switch on Clause Learning. Possible values are: 
		
        * 0 : No Learning
		* 1 : A naive Learning based on the decisions taken
		* 2 : CDCL 
        
        Default value : 0

* -semantic : [binary 1/0] enable Semantic Learning
 
		Default value : 0

* -lazygeneration : [binary 1/0] lazily generate bound literals needed in the learnt clause 
        
        Default value : 0
        

* -orderedexploration : [binary 1/0] Explore the Conflict Graph by following the order of propagation 

        Default value : 0

* -reduce : reduce the Learnt Clause. Possible values {0,1,2}
		
        Default value : 0

* -forgetsize : [positive integer] Used to learn only the nogoods having a size < to this parameter. 
		
        Default value : 0 (no check is performed)

* -forgetbackjump [positive integer] Used to learn only the nogoods leading to a difference of backtrack levels at least equal to this parameter
		
        Default value : 0 (no check is performed)

* -keeplearning [binary 1/0] Used to Enable Learning with Branch and Bound (B&B)
		
        Default value : 0

* -nogoodweight [binary 1/0] Use nogood-based weight within the WeightedDegree Heuristic
		
        Default value : 0


* -simulaterestart {0,1,2} : Used to simulate restart whenever a new bound is found with B&B. Value 2 will initialise the heuristic weights.

        Default value : 0 (no check is performed)
        
        
* -weighthistory : Use the weight history between dichotomy steps and/or B&B. 
 
        * Value  =0 : No history used
        * Value  =1 : The history is updated for each SAT result in the dichotomy but used only to update the weights before starting B&B
        * Value  =2 : The history is updated for each SAT result in the dichotomy but used only to update the weights before each dichotomy step
        * Value  =3 : The history is updated for each SAT result in the dichotomy and used to update the weights before each dichotomy step and before starting B&B
		
        Default value : 0 (no checking is performed)

* -restart : restart policy used with dichotomy : Possible values are :
		
        * "geom" for Geometric Restart 
		* "luby" for Luby Restart 
		
        Default value : "geom"

* -bandbrestart : restart policy used with B&B : Possible values are :
		
        * "geom" for Geometric Restart 
    	* "luby" for Luby Restart 
		
        Default value : "geom"


* -fixedForget : [integer] a fixed number of failure to reach before calling the clause deletion routine.
		
        Default value : 10000

* -fixedlimitSize : [integer] used to trigger clause deletion. If the size of the clause database is at least equal to this parameter then we start forgetting clauses. 
		
        Default value : 2000

* -fixedLearntSize : [integer] A fixed size of the clause database after the reduction.
		
        Default value : 2000


* -maxnogoodsize : [integer] Clauses that have a size at most equal to this parameters will be kept when performing forget. 


* -probforget : [positive integer bounded by 100] The probability that a large clause (w.r.t maxnogoodsize) will be forgotten. By default equal to 100. 



* -forgetfulness : [double] the % of forgetfulness of clauses.

        Default value : 0.0


* -bjmforgetfulness : [double] the % of backjump level forgetfulness of clauses. Very similar to -forgetbackjump, however, with percentage. 
		
        Default value : 0.0 (i.e. no check is performed)

[Mistral-2.0]:https://github.com/ehebrard/Mistral-2.0/
[Taillard instances]:http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/ordonnancement.html 
