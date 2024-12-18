# MASharpe
A Multiobjective Memetic Algorithm for the Markowitz Model based on Informed Decisions

========================================
PISA  (https://sop.tik.ee.ethz.ch/pisa/)
========================================
PISA - A Platform and Programming Language Independent Interface for Search Algorithms.

Permission to use, copy, modify, and distribute this software and its documentation for any purpose, 
without fee, and without written agreement is hereby granted, provided that the above copyright notice and 
the following two paragraphs appear in all copies of this software.

IN NO EVENT SHALL THE SWISS FEDERAL INSTITUTE OF TECHNOLOGY, 
COMPUTER ENGINEERING AND NETWORKS LABORATORY BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, 
INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, 
EVEN IF THE SWISS FEDERAL INSTITUTE OF TECHNOLOGY, COMPUTER ENGINEERING AND NETWORKS LABORATORY HAS 
BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

THE SWISS FEDERAL INSTITUTE OF TECHNOLOGY, COMPUTER ENGINEERING AND NETWORKS LABORATORY, SPECIFICALLY 
DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND 
FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE SWISS 
FEDERAL INSTITUTE OF TECHNOLOGY, COMPUTER ENGINEERING AND NETWORKS LABORATORY HAS NO OBLIGATION TO PROVIDE 
MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 

===============================

MASHARPE - Investment portfolio optimization Markowitz model using Sharpe as a knowledge base

Variator Implementation with PISALib
 
Documentation
  
file: masharpe_documentation.txt
author: Feijoo Colomine D, Carlos Cotta, Antonio Fernandez
last change: 
============================

The Problem
===========

This module is based on the following works and source codes:

-The problems MARKOWITZ, SHARPE described in [Markowitz1952,Sharpe1971].

- The problems DTLZ1 to DTLZ7 described in [DTLZ2002a] (except for
DTLZ6) and [DTLZ2005a] as well as the problem COMET described in the
appendix of [DTLZ2005a]. The numbering in this module is according to
[DTLZ2005a].

- The six two-objective problems proposed in [ZDT2000a],
here named ZDT1 to ZDT6,

- The two-objective problem SPHERE (multi-objective sphere
model), KUR (Kursawe's function) and QV (test problem
proposed by Quagliarella and Vicini), as described in
[ZLT2001a].

Note that the ZDT3 problem was slighlty modified compared to
the original version: in order to conform to the PISA
convention of non-negative objective values, the constant 1
was added to the second function f_2 of ZDT3.

InProceedings{DTLZ2002a,
  author =       {K. Deb and L. Thiele and M. Laumanns
                  and E. Zitzler}, 
  title =        {Scalable Multi-Objective Optimization Test Problems},
  booktitle =    {Congress on Evolutionary Computation (CEC)},
  year =         2002,
  pages     = {825--830},
  publisher = {IEEE Press},
}

@InCollection{DTLZ2005a,
  author =       {K. Deb and L. Thiele and M. Laumanns
                  and E. Zitzler},
  title =        {Scalable Test Problems for Evolutionary
                  Multi-Objective Optimization},
  booktitle =    {Evolutionary Multiobjective Optimization:
                  Theoretical Advances and Applications}, 
  publisher =    {Springer},
  year =         2005,
  isbn =	 {1-85233-787-7},
  chapter =	 6,
  pages =        105--145,
  editor =       {A. Abraham and R. Jain and R. Goldberg},
}

@Article{ZDT2000a,
  author = "E. Zitzler and K. Deb and L. Thiele",
  title = "Comparison of Multiobjective Evolutionary
           Algorithms: Empirical Results", 
  journal = "Evolutionary Computation",
  year = 2000,
  volume =       8,
  number =       2,
  pages =        "173--195"
}

@TechReport{ZLT2001a,
  author = 	 "E. Zitzler and M. Laumanns and L. Thiele",
  title = 	 "{SPEA2}: Improving the {S}trength {P}areto
                  {E}volutionary {A}lgorithm", 
  institution =  "Computer Engineering and Networks Laboratory (TIK), Swiss
		  Federal Institute of Technology (ETH) Zurich",
  year = 	 2001,
  number =	 103,
  address =	 "Gloriastrasse 35, CH-8092 Zurich, Switzerland",
  month =	 "May"
}

@Book{Deb2001a,
  author={K. Deb},
  title={Multi-objective optimization using evolutionary algorithms},
  publisher={Wiley},
  address =  "Chichester, UK",
  year={2001} 
}

@ARTICLE{Markowitz1952,
  author = {Markowitz, Harry M.},
  title = {Portfolio Selection},
  journal = {Journal of Finance},
  year = {1952},
  volume = {7},
  pages = {77-91},
  doi = {10.2307/2975974},
  owner = {MI PC},
  timestamp = {2006.08.03}
}

@Article{Sharpe1971,
  author     = {Sharpe, William F.},
  title      = {Mean-Absolute-Deviation Characteristic Lines for Securities and Portfolios},
  journal    = {Manage. Sci.},
  year       = {1971},
  volume     = {18},
  number     = {2},
  pages      = {B-1–B-13},
  month      = {oct},
  issn       = {0025-1909},
  address    = {Linthicum, MD, USA},
  doi        = {10.1287/mnsc.18.2.B1},
  issue_date = {October 1971},
  numpages   = {14},
  publisher  = {INFORMS},
  url        = {https://doi.org/10.1287/mnsc.18.2.B1},
}


The Variation Operators
=======================

The individuals are represented as real vectors, normalized to the
interval [0,1]. The variation operators are the SBX (simulated binary
crossover) and the polynomial mutation operator (see [Deb2001a]). If
the parameter <variable_recombination_probability> is set to 0 and
<variable_swap_probability> set to 0.5, the recombination is
equivalent to uniform crossover.



The Parameters
==============

This module uses the following values for the common parameters:

alpha   (size of the initial population)
mu      (number of parent individuals)
lambda  (number of offspring individuals, has to be equal to mu)
dim     (number of objectives)

'PISA_cfg' is a sample PISA_configuration file.

MASHARPE takes local parameters which are given in a parameter file. 
The name of this parameter file is passed to the MASHARPE program as
command line argument. 

problem (test problem name, e.g., DTLZ2 or SPHERE OR COLOMBIA20)

seed (seed for random number generator)

number_decision_variables (number of decision variables)

maxgen (Number of objective functions)

outputfile (outputfile)

individual_mutation_probability (probability that a certain individual
undergoes mutation)

individual_recombination_probability (probability that a certain pair
of individuals undergoes recombination) 

variable_mutation_probability (probability that a certain variable in
a given individual is mutated)

variable_swap_probability (probability that a certain pair of
variables is swapped during recombination)

variable_recombination_probability (probability that the SBX
recombination operator is used for a given pair of variables; this
decision is independent from variable_swap_probability)

eta_mutation 20 (distribution index for mutation operator)

eta_recombination 15 (distribution index for recombination operator)

use_symmetric_recombination 1 (switch on or off the symmetry
constraint for recombination)

Parameters inherent to the use of hybrid capabilities:

K: Best Sharpe Individuals Vector Size

random: (seed for random number generator)

numberOfRuns: number of executions of the masharpe algorithm

P_ls: probability of action of local search (P_Ls)

start_generation_number: initial generation for local search or elite memory to work

final_generation_number: final generation for local search or elite memory to work

P_em: Elite memory usage probability using Sharpe (P_em)

max_num_iter_l_s:  maximum number of reviews for local search

number_of_evaluations_max: maximum number of evaluations of the algorithm

initial_population_checker: use the initial population checker (0: Do not use; maximum value 1 and works according to a random function) 

individual_search: local search probability if the initial population checker is activated (corrector_pob_inicial>0)

diversity: Elite population size (diversity<50)

'masharpe_param.txt' is a PISA_parameter file.


Source Files
============

The source code for masharpe is divided into six files.

Four generic files are taken from PISALib:

'variator.{h,c}' is a taken from PISALib. It contains the main
function and all functions implementing the control flow.

'variator_internal.{h,c}' is taken from PISALib. It contains functions
that are called by the functions in the 'variator' part and do the
work in the background (file access etc.). 
  
'variator_user.{h,c}' defines and implements the DTLZ specific
operations.

Additionally, a Makefile, a 'PISA_cfg' file with common parameters and
a 'dtlz_param.txt' file with local parameters used by DTLZ are
contained in the tar file.

For compiling on Windows and Unix (any OS having <unistd.h>) uncomment
the according '#define' in the 'variator.h' file.


Usage
=====

Call MASHARPE with the following arguments:
You can run the defined batch file runhybrid.bat in windows environment (64 bit)

masharpe paramfile filenamebase poll

paramfile: specifies the name of the file containing the local
parameters (e.g. masharpe_param.txt)

filenamebase: specifies the name (and optionally the directory) of the
communication files. The filenames of the communication files and the
configuration file are built by appending 'sta', 'var', 'sel','ini',
'arc' and 'cfg' to the filenamebase. This gives the following names for
the '../PISA_' filenamebase:

../PISA_cfg - configuration file
../PISA_ini - initial population
../PISA_sel - individuals selected for variation
../PISA_var - variated individuals (offspring)
../PISA_arc - individuals in the archive

Caution: the filenamebase must be consistent with the name of
the configuration file and the filenamebase specified for the selector
module.

poll: gives the value for the polling time in seconds (e.g. 0.5).


Output
======

masharpe writes the content of the archive in the last generation to a
specified output file. One individual is written per line using the
following format:

ID (objective 1) (objective 2) ... (objective dim) (decision variable
1) (decision variable 2) ... (decision variable <number_decision_variables>) 



Limitations
===========

This masharpe module can only handle mu == lambda. If an odd number is
chosen for mu and lambda, the last individual in the mating pool can
only undergo mutation, as it has no recombination partner.



Stopping and Resetting
======================

The behaviour in state 7 and 11 is not determined by the interface but
by each variator module specifically. masharpe behaves as follows:

state 7 (= selector terminated): set state to 4 (terminate as well).
state 11 (= selector resetted): set state to 0 (start again).

The user can change the state variable in the sta file using a text
editor, e.g., for stopping both processes or for resetting. MASHARPE 
assumes that the variator is resetted before the selector, i.e., state
8 is present before state 10.
