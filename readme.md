This project is the Dynamic Controllability Checking algorithm for Controllable Conditional Temporal Problem with Uncertainty.
A paper about this work has been published on ICAPS 2017.


# compile
   Need Gurobi and Gurobi license

# run
   cds_stnu is used for calculating
   
   cds_test is used for generating test cases

## TEST Max_delay (MD)
   cds_stnu [filename]

## TEST Dynamic_Controllability (DC Checking)
   cds_stnu -o 5 [filename]
   
---------------------
|   parameters     | DC Checking           | Optimisation (MD)  |
| ------------- |-------------|-----|
| STPU      | -v 1 -o 5 | -v 1 [-o 0] |
| CCTPU     | [-v 0] [-o 0] |   TODO |


# example
   example1.pos is a POS, which is used for testing MD
   
   example2.cst is a CCTPU, which is used for testing DC




