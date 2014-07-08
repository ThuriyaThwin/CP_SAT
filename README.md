[Work in progress]

Description
-----------

This repository is a CDCL version of [Mistral-2.0] dedicated to resolve (Disjunctive) Scheduling Problems. So far, only the standard Job-Shop problem is being (experimentally) used with [Taillard instances]. 

Installation
--------------
To use the Scheduler module, first compile the code using 

```sh
 make scheduler 
```


That's it! Now you can try it with 

```sh
>bin/scheduler data/scheduling/jsp/15x15txt/1.txt 
```
The folder data/scheduling/jsp/15x15txt contains all the 15x15 instances. 

Parameters
----------


Several parameters can be used. These are probably the most important ones used with Clause Learning (I'll update the list later) : 

* -fdlearning : Used to switch on Clause Learning 
	      Values : 
		0 : No Learning
		1 : A naive Learning based on the decisions taken
		2 : CDCL 
		Default Value : 0

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

* -forgetbackjump [positive integer] Used to learn only the nogoods leading to a backtrack level at least equal to this parameter
		Default value : 0 (no check is performed)

* -keeplearning [binary 1/0] Used to Enable Learning with Branch and Bound (B&B)
		Default value : 0

* -nogoodweight [binary 1/0] Use nogood-based weight within the WeightedDegree Heuristic
		Default value : 0
* -simulaterestart {0,1,2} : Used to simulate restart whenever a new bound is found with B&B. Value 2 will initialise the heuristic weights.
		Default value : 0 (no check is performed)
* -weighthistory : Use the weight history between dichotomy steps and/or B&B. 
        Value  =0 : No history used
        Value  =1 : The history is updated for each SAT result in the dichotomy but used only to update the weights before starting B&B
        Value  =2 : The history is updated for each SAT result in the dichotomy but used only to update the weights before each dichotomy step
        Value  =3 : The history is updated for each SAT result in the dichotomy and used to update the weights before each dichotomy step and before starting B&B
		Default value : 0 (no checking is performed)

* -restart : restart policy used with dichotomy : Possible values are :
		"geom" for Geometric Restart 
		"luby" for Luby Restart 
		Default value : "geom"

* -bandbrestart : restart policy used with B&B : Possible values are :
		"geom" for Geometric Restart 
		"luby" for Luby Restart 
		Default value : "geom"

* -forgetfulness : [double] the % of forgetfulness of clauses.
		Default value : 0.0
* -bjmforgetfulness : [double] the % of backjump level forgetfulness of clauses. Very similar to -forgetbackjump, however, with pourcentage. 
		Default value : 0.0 (i.e. no check is performed)

[Mistral-2.0]:https://github.com/ehebrard/Mistral-2.0/
[Taillard instances]:http://mistic.heig-vd.ch/taillard/problemes.dir/ordonnancement.dir/ordonnancement.html 
