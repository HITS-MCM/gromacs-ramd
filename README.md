[![Build Status](https://jenkins.h-its.org/buildStatus/icon?job=MCM/HITS-MCM/gromacs-ramd/main)](https://jenkins.h-its.org/job/MCM/job/HITS-MCM/job/gromacs-ramd/job/main/)

# Random Acceleration Molecular Dynamics (RAMD)

## Background

Random Acceleration Molecular Dynamics (RAMD) is a method to carry out
molecular dynamics simulations with an additional randomly oriented force
applied to a molecule in the system.

## Installation

See http://manual.gromacs.org/documentation/current/install-guide/index.html

## HPC deployment with EasyBuild

The [EasyBuild](http://easybuilders.github.io/easybuild/) recipe is available [here](https://github.com/BerndDoser/easybuild-easyconfigs/blob/hits/easybuild/easyconfigs/g/GROMACS/GROMACS-2020.5-RAMD-2.0-fosscuda-2019b.eb).

## Usage

Please use following mdp options:

 * ramd

   RAMD will be applied.

 * ramd-seed

   Seed for random direction generator

 * ramd-eval-freq

   This parameter affect absolute dissociation time but have less
   effect on the relative dissociation times of different compounds. It is
   recommended to use default value.

 * ramd-force-out-freq

   This ramd parameter resets pull-nstxout and pull-nstfout.

 * ramd-ngroups

   The number of ramd groups defining the ligand-receptor pair.
   Below only the pull options for group 1 are given,
   further groups simply increase the group index number.

 * ramd-group1-force

   The force constant in kJ/mol/nm. Default value is 600 kJ/mol/nm.
   For a set of compounds with the dissociation rate expected to vary
   within the range of 0.1-0.0001 1/s, a random force magnitude of 600
   kJ/mol/nm can be applied. If necessary, the force magnitude can be
   adjusted according to the longest and shortest dissociation time
   observed in simulations. The upper threshold of the force magnitude is
   determined by the fast-dissociated compounds, whose dissociation time
   should be longer than 100 ps. The lower threshold of the force magnitude
   depends on the computation facilities available.

 * ramd-group1-r-min-dist

   This parameter affect absolute dissociation time but have less
   effect on the relative dissociation times of different compounds. It is
   recommended to use default value.

 * ramd-group1-max-dist

   This value has to be adjusted for the system studied: no
   protein-ligand contacts should be observed in the last snapshot of a
   dissociation trajectory. Usually 4 nm is enough, but in the case
   of a long dissociation channel (as in many membrane proteins) maxDist must be
   increased accordingly. Method performance is not very sensitive to the
   upper limit of this parameter since motion of the free ligand due to the
   external force is very fast (i.e. the last part of the trajectory, where the
   ligand does not interact with the protein, usually has a negligible contribution
   to the observed dissociation time).

 * ramd-group1-receptor

   Receptor for the first RAMD group. Default name is 'Protein'.

 * ramd-group1-ligand

   Ligand for the first RAMD group. Default name is 'INH'.

 * ramd-group1-receptor-pbcatom

   The value will be forwarded to the associated pull group of the receptor.
   Default value is 0, which takes the middle atom (number wise).

 * ramd-group1-ligand-pbcatom

   The value will be forwarded to the associated pull group of the ligand.
   Default value is 0, which takes the middle atom (number wise).

 * ramd-pbc-ref-prev-step-com

   The value will be forwarded to pull-pbc-ref-prev-step-com. Default value is 'yes'.

 * ramd-connected-ligands

   If ‘yes’, the trajectory will be terminated when all ligands have left the radius.
   If one ligand leaves the radius, its last assigned force will continue to be
   applied until the simulation end or radius re-entry.
   If ‘no’, this should revert to standard RAMD for multiple disconnected ligands.
   Each ligand is subject to a RAMD force until the individual ligand has left the
   dissociation radius. The simulation stops when all ligands have left the
   dissociation radius. Default value is 'yes'.
