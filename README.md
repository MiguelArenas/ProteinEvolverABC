# ProteinEvolverABC
Estimation of Recombination and Substitution rates in alignments of protein sequences by approximate Bayesian computation

The package ProteinEvolverABC is a computer framework to estimate recombination and substitution rates in multiple alignments of protein sequences by approximate Bayesian computation. The framework is based on a special version of the simulator ProteinEvolver that implements the coalescent with recombination (including a variety of migration models, demographics and user-specified populations/species trees) and protein evolution under diverse substitution models. A total of 16 summary statistics can be computed from the alignments collecting information about the evolutionary process. Then, ProteinEvolverABC can estimate the parameters under approximate Bayesian computation (ABC) through the rejection and regression methods. 
The framework has been validated with extensive computer simulations. 

The package is implemented in C, perl and R and can run on Linux and Mac OS. It includes a graphical user interface written in Java. Conveniently the simulations can, optionally, run in parallel on a user-specified number of processors. 

The package includes detailed documentation (which is advised to read before using the program) and practical examples.


Citation:

Arenas M. 2022. ProteinEvolverABC: Coestimation of Recombination and Substitution Rates in Protein Sequences by approximate Bayesian computation. Bioinformatics, 38(1), 58–64.


Miguel Arenas. 2021-2022. miguelmmmab@gmail.com - marenas@uvigo.es
