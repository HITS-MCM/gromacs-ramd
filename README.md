[![Build Status](https://jenkins.h-its.org/buildStatus/icon?job=MCM/HITS-MCM/gromacs-ramd/release-2019)](https://jenkins.h-its.org/job/MCM/job/HITS-MCM/job/gromacs-ramd/job/release-2019/)

# Random accelerated Molecular Dynamics (RAMD)

## Background

Random Acceleration Molecular Dynamics (RAMD) is a method to carry out
molecular dynamics simulations with an additional randomly oriented force
applied to a molecule in the system.

## Installation

See http://manual.gromacs.org/documentation/2019/install-guide/index.html

## Usage

Please use following mdp options:

 * ramd: Switch on RAMD
 * ramd-receptor
 * ramd-ligand
 * ramd-force: The force constant

   For a set of compounds with the dissociation rate expected to vary
   within the range of 0.1-0.0001 1/s, a random force magnitude of 600
   kJ/mol/nm can be applied. If necessary, the force magnitude can be
   adjusted according to the longest and shortest dissociation time
   observed in simulations. The upper threshold of the force magnitude is
   determined by the fast-dissociated compounds, whose dissociation time
   should be longer than 100 ps. The lower threshold of the force magnitude
   depends on the computation facilities available.

 * ramd-seed

   Seed for random direction

 * ramd-eval-freq

   This parameter affect absolute dissociation time but have less
   effect on the relative dissociation times of different compounds. It is
   recommended to use default value.

 * ramd-r-min-dist

   This parameter affect absolute dissociation time but have less
   effect on the relative dissociation times of different compounds. It is
   recommended to use default value.

 * ramd-force-out-freq

   This ramd parameter resets pull-nstxout and pull-nstfout.

 * ramd-max-dist

   This value has to be adjusted for the system studied: no
   protein-ligand contacts should be observed in the last snapshot of a
   dissociation trajectory. Usually 4 nm is enough, but in the case
   of a long dissociation channel (as in membrane proteins) maxDist must be
   increased accordingly. Method performance is not very sensitive to the
   upper limit of this parameter since motion of the free ligand due to the
   external force is very fast (i.e. the last part of the trajectory, where
   ligand does not interact with the protein, has a negligible contribution
   to the observed dissociation time.

