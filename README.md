# ProteinEvolverABC2
Estimation of Recombination and Substitution rates in alignments of protein sequences by approximate Bayesian computation

The package ProteinEvolverABC is a computer framework to estimate recombination and substitution rates in multiple alignments of protein sequences by approximate Bayesian computation. The framework is based on a special version of the simulator ProteinEvolver that implements the coalescent with recombination (including a variety of migration models, demographics and user-specified populations/species trees) and protein evolution under diverse substitution models. A total of 20 summary statistics can be computed from the alignments collecting information about the evolutionary process. Then, ProteinEvolverABC can estimate the parameters under approximate Bayesian computation (ABC) through the rejection and regression methods.
The second version of the framework allows for the estimation of these parameters accounting for evolutionary constraints from the protein structure.
The framework has been validated with extensive computer simulations. 

The package is implemented in C, perl and R and can run on Linux and Mac OS. It includes a graphical user interface written in Java. Conveniently the simulations can, optionally, run in parallel on a user-specified number of processors. 

The package includes detailed documentation (which is advised to read before using the program) and practical examples.

o	The first version of this framework (now available from Releases, https://github.com/MiguelArenas/ProteinEvolverABC/releases) was supported by the Spanish Ministry of Economy and Competitiveness and Ministry of Science and Innovation through the Grants [RYC- 2015-18241] and [PID2019-107931GA-I00/AEI/10.13039/501100011033]. Arenas M. 2022. ProteinEvolverABC: Coestimation of Recombination and Substitution Rates in Protein Sequences by approximate Bayesian computation. Bioinformatics, 38(1), 58–64.

o	The second (current) version of this framework was partially supported by the grants [PID2023-151032 NB-C21] and [PID2023-151032NB-C22] funded by the Spanish Ministry of Science MCIU/the Spanish Agency of Research AEI/10.13039/501100011033/ and FEDER, UE. This version is also available from Releases, https://github.com/MiguelArenas/ProteinEvolverABC/releases


Elena Pazos-Linares, Cristina Landa, Ugo Bastolla, Miguel Arenas
2020-2026. marenas@uvigo.es

