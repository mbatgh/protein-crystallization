# protein-crystallization
A fast solver for population balance equations, that can be used to model the impact of process parameter settings on the results of anti-solvent crystallization.

## Introduction
This is a code I developed to model and understand the impact of process parameters on the yield of anti-solvent crystallization from aqueous solution. As it is rather fast, it could be used as part of PAT equipment to control crystallization as purification step during the manufacture of proteins. Crystal nucleation, growth, and aggregation are accounted for quantitatively via a population balance equation solver with appropriate kernels. The solver is based on functions provided by the boost ODEINT library. Material parameters, such as nucleation, growth, and aggregation rates need to be provided as input. These parameters can be established using lab scale experiments.

## Compilation
The program is written in vanilla C++. So far it has only been tested on Linux (Debian and Ubuntu) with gcc. However, it should work on most reasonably recent Linux distros, as long as boost is installed. Compilation is straight forward. Clone the repository, cd to the downloaded folder, and enter:

make

## Usage
Compilation provides a binary called aspbe. Copy it to a folder in your PATH, and executed on the command line, as in:

aspbe -h

This provides a list of input parameters. To run a simulation you need to provide at least a configuration file (-f), and a run ID (-i). THe former is a text file containing parameters such as initial seed, growth rate, etc. The latter is a string used for naming of output files. Examples for configuration files are in the test folder.
Detailed documentation is not available at this point in time.
