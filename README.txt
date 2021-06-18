Documentation of the module main.py used for the experimentation of: Towards a 
Pareto Front Shape Invariant Multi-Objective Evolutionary Algorithm Using 
Pair-Potential Functions.

NAME
       main.py - test a multi-objective evolutionary algorithm (MOEA)

SYNOPSIS
       main.py ARG1 ARG2 ARG3 ARG4 ARG5 ARG6 ARG7
       main.py OPTION

DESCRIPTION
       This module is used to test a MOEA on a selected multi-objective problem 
       (MOP) with a given number of objective functions for a specific number 
       of independent runs. The required arguments are described as follows:

       ARG1
              It must be a valid MOEA name. The valid MOEA names are: NSGA-III, 
              NSGA-III-RSE, NSGA-III-GAE, NSGA-III-COU, NSGA-III-PT, 
              NSGA-III-MPT, NSGA-III-GPT, and NSGA-III-KRA.

       ARG2
              It must be an integer greater than zero. It represents the number 
              of divisions per objective function for the boundary layer used
              for the generation of a weight vector-based reference set.

       ARG3
              It must be an integer greater than or equal to zero. It 
              represents the number of divisions per objective function for the 
              inner layer used for the generation of a weight vector-based 
              reference set. The inner layer is not generated when a value of 
              zero is given.

       ARG4
              It must be a valid MOP name. The valid MOP names are: DTLZ1, 
              DTLZ2, DTLZ3, DTLZ4, DTLZ5, DTLZ6, DTLZ7, DTLZ1_MINUS, 
              DTLZ2_MINUS, DTLZ3_MINUS, DTLZ4_MINUS, DTLZ5_MINUS, DTLZ6_MINUS, 
              and DTLZ7_MINUS.

       ARG5
              It must be an integer greater than one. It represents the number 
              of objective functions of the MOP.

       ARG6
              It must be an integer greater than or equal to zero. It 
              represents the maximum number of generations for the execution of
              the MOEA.

       ARG7
              It must be an integer greater than zero. It represents the number 
              of independent runs that the MOEA will be tested on the selected 
              MOP.

       The following option can be used:

       --help 
              Display this help and exit.

REQUIREMENTS
       A computer with the installation of Python 3.8 is needed. The modules 
       numpy, matplotlib, and scipy are required.

EXAMPLE
       For running the module main.py, go to NSGA-III-K/ and write:

       IPython console users:
              %run main.py NSGA-III-RSE 12 0 DTLZ1_MINUS 3 400 1

       Windows users:
              python main.py NSGA-III-RSE 12 0 DTLZ1_MINUS 3 400 1

       Linux users:
              python3 main.py NSGA-III-RSE 12 0 DTLZ1_MINUS 3 400 1

       The previous line executes the module main.py to test the NSGA-III-RSE 
       on the DTLZ1_MINUS with three objective functions for one independent 
       run. A boundary layer with twelve divisions per objective function is
       used for the weight vector-based reference set and the inner layer is 
       not required. The maximum number of generations is set to four hundred 
       generations.

RESULTS
       On success, the output files containing the approximation sets are 
       generated in NSGA-III-K/Results/.
