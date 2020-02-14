#***********************************************************************************
#                                                                                  *
# Random Acceleration Molecular Dynamics (RAMD)                                    *
# Implementation for NAMD v2.12                                                    *
# September 2017                                                                   *
#                                                                                  *
# Copyright (c) 2017, EML Research gGmbH, Heidelberg, Germany                      * 
# Authors: Vlad Cojocaru, Stefan Richter, Daria Kokh, Rebecca Wade                 *
# Email: mcmsoft@h-its.org                                    *
#                                                                                  *
# The first Tcl script to run RAMD in NAMD (v2.5+) was written by Harish Vashisth  *
# Ref: Vashisth H et al, Biophys J. 2008 Nov 1;95(9):4193-204. Epub 2008 Aug 1     *
#                                                                                  *
# The Vashisth script inspired some lines of this script, but mainly this script   *
#    is based on the fortran code written for AMBER 8                              *
#                                                                                  *
# The structure of this script is inspired by the Adaptive Biasing Force module    *
#    distributed with NAMD 2.6                                                     *
#                                                                                  *
# The original RAMD method is described in:                                        *
#    Ref1: Luedemann,S.K.,Lounnas,V.and R.C.Wade.,                                 *
#          J Mol Biol, 303:797-811 (2000)                                          *
#    Ref2: Schleinkofer,K.,Sudarko,Winn,P.,Luedemann,S.K.and R.C.Wade,             *
#          EMBO Reports, 6, 584-589 (2005)                                         *
#                                                                                  *
# The present script is a slightly changed version of ramd-5_script.tcl            *
#  adjusted to use a force magntude asan input for Kof simulatios                  *
#                                                                                  *
# Disclaimer:  This script is for research purposes only. EML Research does not    *
#              assume any responsibility for the software or its use.              * 
#                                                                                  *
# !!! Quantitative reproducibility of the results obtained with AMBER 8            *
#     is not possible due to a numerical errors and such like                      * 
#                                                                                  *
# !!! This script is under development.                                            *
#                                                                                  *
#   The script along with usage examples is available at                           *
#   https://www.h-its.org/en/research/mcm/software/                                *
#***********************************************************************************
namespace eval ::RAMD {
# set ramdfilename "ramd.log"
 set ramdfileid [open $ramdfilename "w"]
 fconfigure $ramdfileid -blocking 0
 if {$namdVersion > 2.10} {
   enabletotalforces
 }
 puts $ramdfileid "RAMD:"
 puts $ramdfileid "RAMD:   -------------------------------------------------------------------  "
 puts $ramdfileid "RAMD:   Random Acceleration Molecular Dynamics Simulation version $version"
 puts $ramdfileid "RAMD:   -------------------------------------------------------------------  "
 puts $ramdfileid "RAMD:"
 #***** Assign default values for the parameters not specified in the configuration file
 foreach option [array names defaults] {
  if {! [info exists $option]} {
   set $option $defaults($option)
   puts $ramdfileid [format "RAMD: %25s = %s" $option [expr $$option]]
  } elseif { [info exists $option] } {
   puts $ramdfileid [format "RAMD: %25s = %s" $option [expr $$option]]
  }
 }
 #***** Check if mandatory parameters are specified in the configuration file
 foreach var $mandatory {
  if {! [info exists $var]} {
   error "RAMD: Mandatory parameter '$var' is not set -- cannot start RAMD"
  } else {
   puts $ramdfileid [format "RAMD: %25s = %s" $var [expr $$var]]
  }
 }
 #***** Check if 'forceOutFreq' is equal to 1; exit with error if that's the case
 if { $forceOutFreq == 1 } { error "RAMD: ERROR: 'forceOutFreq' parameter may not be 1" } 

 #***** Check if 'mdSteps' is specified in the configuration file

 #***** Performed pure RAMD if 'mdSteps' = 0
 if { $mdSteps == 0 } { 
  
  #***** Check that the number of ramd steps is a multiple of 'forceOutFreq'; exit with error if not
  set r [expr "$ramdSteps % $forceOutFreq"]
  if { $r != 0 } { error "RAMD: ERROR: The number of RAMD steps is not a multiple of 'forceOutFreq'" } 
 
  puts $ramdfileid "RAMD: Pure RAMD simulation is performed" 
  
  #***** If 'mdSteps' is 0 and "mdStart" is yes, give a warning
  if { $mdStart == "yes" } { 
   puts $ramdfileid "RAMD: WARNING: 'mdStart' has no meaning for pure RAMD simulation; it will be ignored" 
  }

  #***** If 'mdSteps' is 0 and "rMinMd" is set, give a warning
  if { [info exists rMinMd] } {
   puts $ramdfileid "RAMD: WARNING: 'rMinMd' specified while 'mdSteps' is 0"
   puts $ramdfileid "RAMD: WARNING: For combined RAMD-MD simulation 'mdSteps' must be greater than 0"
   puts $ramdfileid "RAMD: WARNING: Ignore 'rMinMd' and perform pure RAMD simulation"
  }
  
 }

 #***** Perform combined RAMD with MD simulation if 'mdSteps' is not 0 and 'rMinMd' is specified
 if { $mdSteps != 0  } { 
  error "To run ramd setting mdSteps is not supported"

 }
     
 #***** Make a list of all the atoms on which the force will be applied
 set ramdAtoms {}
 for { set i $firstRamdAtom } { $i <= $lastRamdAtom } { incr i } { lappend ramdAtoms $i }
 puts $ramdfileid "RAMD: Atoms subject to the random acceleration are: $ramdAtoms"
 foreach ramdAtom $ramdAtoms { addatom $ramdAtom }
 #***** Define a group of the ligand atoms; the force will be applied on the center of mass of this group
 set ramdGroup [ addgroup $ramdAtoms ]

 #***** Define a group containing all protein atoms
 set protAtoms {}
 for { set i $firstProtAtom } { $i <= $lastProtAtom } { incr i } { lappend protAtoms $i }
 foreach protAtom $protAtoms { addatom $protAtom }
 set protGroup [ addgroup $protAtoms ]

 #***** Some variable initialization 
 set timeStep 0; set exitFlag 0; 
 set prevLigCOM "0.0 0.0 0.0"; set prevProtCOM "0.0 0.0 0.0"; set prevDist 0;

 #***** Initialization of simulation flags
 if { $mdSteps == 0 } {
  set ramdFlag 1; set mdFlag 0; set ramdStep 0; set mdStep 0;
 }
 
} ;# namespace
 
#***** In root namespace (::) for all procedures we have to add the following procedure definition
#proc veclength {v} {
# return [expr {sqrt([veclength2 $v])}]
#}
#***** Source the vectors and matrices procedures from VMD
#source $RAMD::RAMDdir/vectors.tcl

#***********************************************************
# PROCEDURE TO GENERATE RANDOMLY ORIENTED ACCELERATION 
#***********************************************************
proc genRandAccel { timeStep } {
namespace eval ::RAMD {

 set pi 3.141592653589793238462643383
# set pi [expr "2.0*asin(1.0)"]

 #***** Generate new random orientation of the ramd force
# set randTheta [expr "rand()"]
# set randPsi [expr "rand()"]
 set theta [expr "2*$pi*rand()"]
 set psi [expr "$pi*rand()"]
 set sinpsi [expr "sin($psi)"]
 set rx [expr "cos($theta)*$sinpsi"]
 set ry [expr "sin($theta)*$sinpsi"]
 set rz [expr "cos($psi)"]
 set r "$rx $ry $rz"
 set lenr [veclength $r]

 # Acceleration is given in kcal/mol*A*amu  in the NAMD configuration file (multiply with 418.68 to get A/ps^2)
 set vecAccel [vecscale [expr "$forceRAMD"] $r ]
  
 return 
 
} ;# namespace
} ;# proc genRandAccel {timestep}


#*****************************************************************************
# PROCEDURE TO EVALUATE THE DISTANCE TRAVELLED BY THE LIGAND IN N RAMD STEPS
#*****************************************************************************
proc evalWalkDist { timeStep prevLigCOM prevProtCOM currLigCOM currProtCOM } {
namespace eval ::RAMD {
 
 #***** Compute the relative position of the ligand com with regard to the protein com
 set prevRelLigCOM [ vecsub $prevLigCOM $prevProtCOM ]
 set currRelLigCOM [ vecsub $currLigCOM $currProtCOM ]
 
 #***** Compute the distance travelled by the ligand com during a ramd or md stint
 set vecWalkDist [vecsub $currRelLigCOM $prevRelLigCOM]
 set walkDist [veclength $vecWalkDist]

# set vecWalkDistX [lindex $vecWalkDist 0]
# set vecWalkDistY [lindex $vecWalkDist 1]
# set vecWalkDistZ [lindex $vecWalkDist 2]
 
 return  

} ;# namespace
} ;# proc evalWalkDist
  

#**************************************************************
# PROCEDURE TO APPLY THE FORCE WHICH IS CALLED EVERY TIME STEP 
#**************************************************************
proc calcforces {} {
namespace eval ::RAMD {
 #***** Terminate NAMD if the ligand has exited from the protein
 if { $exitFlag == 1 } {
  puts $ramdfileid "EXIT: $timeStep  > MAX DISTANCE LIGAND COM - PROTEIN COM REACHED"
  puts $ramdfileid "EXIT: $timeStep  > LIGAND EXIT EVENT DETECTED: STOP SIMULATION"
  puts $ramdfileid "EXIT: $timeStep  > EXIT NAMD"
  fconfigure $ramdfileid -blocking 1
  close $ramdfileid
  puts "EXIT: $timeStep  > MAX DISTANCE LIGAND COM - PROTEIN COM REACHED"
  puts "EXIT: $timeStep  > LIGAND EXIT EVENT DETECTED: STOP SIMULATION"
  puts "EXIT: $timeStep  > EXIT NAMD"
  set process [pid]
  puts "Killing $process"
  exec kill -15 $process
#  exec kill -15 $$
  exit 0
#  exec /bin/bash qdel_jobid.sh
 } 

 if { [ array exists coords ] } { array unset coords }
# if { [ array exists masses ] } { array unset masses }
# if { [ array exists extForces ] } { array unset extForces }
 if { [ array exists totForces ] } { array unset totForces }

 #***** Load coordinates for all the atoms and groups defined
 loadcoords coords
 #***** Load masses for all the atoms and groups defined
 if  {[array exists masses] == 0} { loadmasses masses }
 #***** Calculate the mass of the ligand
 if {[info exists ligMass] == 0} {
    set ligMass 0
    foreach ramdAtom $ramdAtoms {
      set ligMass [expr $ligMass + $masses($ramdAtom)]
    }
    puts $ramdfileid "RAMD: calculated ligand mass $ligMass"
    foreach ramdAtom $ramdAtoms {
      set relAtomMass($ramdAtom) [expr "$masses($ramdAtom)/$ligMass"]
    }
 }

 #***** Load external forces from previous time step for all the atoms and groups defined
# loadforces extForces
 #***** Load total forces from previous time step for all the atoms and groups defined
 loadtotalforces totForces 
 

 #***** Calculate the position of protein and ligand COM
 set protCOM "$coords($protGroup)"
 set ligCOM "$coords($ramdGroup)"
  
 #***** Initialize ramd simulation or combined ramd-md simulation that begins with ramd
 if { $timeStep == 0 } {
  
  expr "srand($ramdSeed)"
  
  set vMin [ expr "($rMinRamd)/($ramdSteps)" ]
  if { $mdSteps == 0 } { 
   puts $ramdfileid "RAMD: $timeStep  ***** INITIALIZE RAMD SIMULATION *****" 
  } else { 
   puts $ramdfileid "RAMD: $timeStep  ***** INITIALIZE COMBINED RAMD-MD SIMULATION *****" 
  }
  puts $ramdfileid "RAMD: $timeStep     >>> minimum travelled distance (A): $rMinRamd"
  puts $ramdfileid "RAMD: $timeStep     >>> minimum velocity (A/fs): $vMin"

  #***** Initial com positions
  set currLigCOM $ligCOM; set currProtCOM $protCOM  
  puts $ramdfileid "RAMD: $timeStep     >>> LIGAND COM IS: $currLigCOM"
  puts $ramdfileid "RAMD: $timeStep     >>> PROTEIN COM IS: $currProtCOM"

  #***** Evaluate initial distance between ligand com and protein com
  set currDist [veclength [vecsub $currLigCOM $currProtCOM]]
  puts $ramdfileid "RAMD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: DIST = $currDist"

  #***** Generate an initial orientation for the acceleration to be applied when ramd is switched on
  genRandAccel $timeStep
  puts $ramdfileid "RAMD: $timeStep     >>> INITIAL RANDOM DIRECTION: $r :: ||r|| = $lenr"

  puts $ramdfileid "RAMD: $timeStep  ***** START WITH $ramdSteps STEPS OF RAMD SIMULATION *****"
 
  #***** Reset the positions of the ligand and protein COMs and the distance ligand com - protein com
  set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

  incr timeStep
  return  

 } 

 
 #***** Perform ramd simulation 
 if { $timeStep != 0 && $exitFlag == 0  } {

  #***** Count ramd steps
  incr ramdStep
#  if { $debugLevel != 0 } {
#   puts $ramdfileid "RAMD DEBUG: TIMESTEP IS: $timeStep; RAMD STEP IS: $ramdStep; MD STEP IS: $mdStep"
#  }
    
  #***** Define and apply the force to each atom of the ligand
#  foreach ramdAtom $ramdAtoms {
#   set atomMass "$masses($ramdAtom)"
#   set atomForce [ vecscale $atomMass $vecAccel ]
#   set atomForce [ vecscale "$relAtomMass($ramdAtom)" $vecAccel ]
#   addforce $ramdAtom $atomForce
#   if { $debugLevel != 0 } { 
#      set atomForceValue [ veclength "$atomForce" ]
#      puts "RAMD DEBUG: ATOM $ramdAtom: MASS $atomMass: ADD FORCE $atomForceValue" 
#      unset atomForceValue
#   }
#   unset atomForce
#   unset atomMass
#  }

  #***** Define the force vector that is applied to the center of mass of the ligand
#  set force [vecscale $ligMass $vecAccel]
#  set force $vecAccel
#  set fx [lindex $force 0]
#  set fy [lindex $force 1]
#  set fz [lindex $force 2]
#  set newaccel [ vecscale [expr "1.0/$ligMass"] $vecAccel ]
 
  #***** Check the magnitude of the force
#  set f [expr "sqrt((($fx)*($fx)+($fy)*($fy)+($fz)*($fz)))"]
   set f [veclength $vecAccel]
  
  #***** Set flag for writing force output
  if { $forceOutFreq != 0 } { set outputFlag [expr "$ramdStep % $forceOutFreq"] }
  
  #***** Write force output every 'forceOutFreq' steps 
  if {  $outputFlag == 0 } {
      
   #***** Write the force vector and direction
   puts $ramdfileid "RAMD FORCE: $timeStep  > LIGAND COM is: $ligCOM\nRAMD FORCE: $timeStep  > PROTEIN COM IS $protCOM\nEXTERNAL FORCE VECTOR (F): $vecAccel; ||F|| = $f\nRAMD FORCE: $timeStep  > EXTERNAL FORCE DIRECTION (r): $r; ||r|| = $lenr"
 #  puts $ramdfileid "RAMD FORCE: $timeStep  > EXTERNAL ACCELERATION VECTOR (F): ||Accel|| = $newaccel"
   }
   #***** Calculate external and total forces acting on the ligand
#   set totLigForceX 0; set totLigForceY 0; set totLigForceZ 0; set totLigForceV 0;
#   set extLigForceX 0; set extLigForceY 0; set extLigForceZ 0; set extLigForceV 0;
   set totLigForce "0 0 0"
#   set extLigForce "0 0 0"
   foreach ramdAtom $ramdAtoms {
     set atomForce [ vecscale "$relAtomMass($ramdAtom)" $vecAccel ]
     addforce $ramdAtom $atomForce
#    set atomMass "$masses($ramdAtom)"
    set totAtomForce "$totForces($ramdAtom)"
#    set totAtomForceX [ lindex $totAtomForce 0 ]
#    set totAtomForceY [ lindex $totAtomForce 1 ]
#    set totAtomForceZ [ lindex $totAtomForce 2 ]
#    set totAtomForceV [ veclength "$totAtomForce" ]
    set totLigForce [vecadd $totLigForce $totAtomForce]
#    set totLigForceX [ expr "$totLigForceX + $totAtomForceX" ]
#    set totLigForceY [ expr "$totLigForceY + $totAtomForceY" ]
#    set totLigForceZ [ expr "$totLigForceZ + $totAtomForceZ" ]
#    set totLigForceV [ expr "$totLigForceV + $totAtomForceV" ]
#    set extAtomForce "$extForces($ramdAtom)"
#    set extAtomForceX [ lindex "$extAtomForce" 0 ]
#    set extAtomForceY [ lindex "$extAtomForce" 1 ]
#    set extAtomForceZ [ lindex "$extAtomForce" 2 ]
#    set extAtomForceV [ veclength "$extAtomForce" ]
#    set extLigForceX [ expr "$extLigForceX + $extAtomForceX" ]
#    set extLigForceY [ expr "$extLigForceY + $extAtomForceY" ]
#    set extLigForceZ [ expr "$extLigForceZ + $extAtomForceZ" ]
#    set extLigForce [vecadd $extLigForce $extAtomForce]
#    set extLigForceV [ expr "$extLigForceV + $extAtomForceV" ]
#    if { $debugLevel != 0 } { 
#     puts $ramdfileid "RAMD DEBUG: ATOM $ramdAtom: MASS $atomMass: EXT FORCE $extAtomForceV"
#     puts $ramdfileid "RAMD DEBUG: ATOM $ramdAtom: MASS $atomMass: TOT FORCE $totAtomForceV"
#    }    
#    unset atomMass
#    unset totAtomForceX; unset totAtomForceY; unset totAtomForceZ; unset totAtomForceV; unset totAtomForce
#    unset extAtomForceX; unset extAtomForceY; unset extAtomForceZ; unset extAtomForceV; unset extAtomForce
   }
#   set totLigForce "$totLigForceX $totLigForceY $totLigForceZ"
#   set extLigForce "$extLigForceX $extLigForceY $extLigForceZ"
 
   #***** Write external forces acting on the ligand com for debugging purposes
#   if { $debugLevel !=0 } {
#    set extLigForceV [veclength $extLigForce]
#    puts $ramdfileid "RAMD DEBUG: $timeStep  > EXTERNAL FORCE ON THE LIGAND COM IS: $extLigForce ($extLigForceV)"
#    unset extLigForceV; 
#   }

   set totLigForceV [veclength $totLigForce]
   #***** Write total forces acting on the ligand com
  if {  $outputFlag == 0 } {

   puts $ramdfileid "RAMD FORCE: $timeStep  > TOTAL FORCE ON THE LIGAND COM IS: $totLigForce ($totLigForceV)"
  }   
   unset totLigForce; #unset extLigForce; 
   unset totLigForceV;

 
#  elseif { !  [ array exists totForces ] && $outputFlag == 0 } {
  
#   error "RAMD: $timeStep  > ERROR: EXTERNAL FORCES NOT PRESENT DURING RAMD STEP: EXIT NAMD"
  
#  }

  #***** Set flag for evaluating ramd simulation
  set evalRamdFlag [expr "$ramdStep % $ramdSteps"]

  #***** Evaluate ramd stint
  if { $evalRamdFlag == 0 } {
  
   puts $ramdfileid "RAMD: $timeStep  ***** EVALUATE $ramdSteps RAMD STEPS AT TIMESTEP $timeStep *****"

   #***** com positions
   set currLigCOM $ligCOM; set currProtCOM $protCOM
#   if { $debugLevel !=0 } {
#    puts $ramdfileid "RAMD DEBUG: $timeStep  > PREVIOUS LIGAND COM IS: $prevLigCOM"
#    puts $ramdfileid "RAMD DEBUG: $timeStep  > PREVIOUS PROTEIN COM IS: $prevProtCOM"
#    puts $ramdfileid "RAMD DEBUG: $timeStep  > CURRENT LIGAND COM IS: $currLigCOM"
#    puts $ramdfileid "RAMD DEBUG: $timeStep  > CURRENT PROTEIN COM IS: $currProtCOM"
#   }
  
   #***** Evaluate distance between ligand com and protein com
   set currDist [veclength [vecsub $currLigCOM $currProtCOM]]   
   #***** Evaluate the change in the distance between the protein and the ligand com during the ramd stint 
   set diffDist [expr "${currDist}-${prevDist}"]
   puts $ramdfileid "RAMD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: $currDist; IT CHANGED BY $diffDist"

   #***** Check if the ligand has exited the protein
   if { $currDist >= $maxDist } { set exitFlag 1; return }
   
   #***** Compute the distance travelled by the ligand com during the ramd stint
   evalWalkDist $timeStep $prevLigCOM $prevProtCOM $currLigCOM $currProtCOM
 
   #***** Evaluate whether a new force direction will be generated
   if { $walkDist <= $rMinRamd } {    

    genRandAccel $timeStep

    puts $ramdfileid "RAMD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (< $rMinRamd)\nRAMD: $timeStep     >>> CONTINUE WITH $ramdSteps STEPS OF RAMD SIMULATION\nRAMD: $timeStep     >>> CHANGE ACCELERATION DIRECTION TO: $r; ||r|| = $lenr"
    
    #***** Reset the ramd step count
    #***** Reset the positions of the ligand and protein COMs
    set ramdStep 0; set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"
 
    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
  
   } elseif { $walkDist > $rMinRamd } {
    puts $ramdfileid "RAMD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (> $rMinRamd)\nRAMD: $timeStep     >>> CONTINUE WITH $ramdSteps STEPS OF RAMD SIMULATION\nRAMD: $timeStep     >>> KEEP PREVIOUS ACCELERATION DIRECTION: $r; ||r|| = $lenr"
   
    #***** Reset the positions of the ligand and protein COMs
    set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"
 
    #***** Increment time step and go to the next time step right now
    incr timeStep
    return
  
   }  
   #***** Ensure that the positions of the ligand and protein COMs are reset after the evaluation
   set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

   #***** Increment time step and go to the next time step right now
   incr timeStep
   return
  
  }

  #***** Unset the force
#  if { [info exists force] } { unset force }
#  if { [info exists fx] } { unset fx }
#  if { [info exists fy] } { unset fy }
#  if { [info exists fz] } { unset fz }
#  if { [info exists f] } { unset f }
#  if { [info exists totLigForce] } { unset totLigForce }
#  if { [info exists extLigForce] } { unset extLigForce }
#  if { [info exists totLigForceV] } { unset totLigForceV }
#  if { [info exists extLigForceV] } { unset extLigForceV }
  
  #***** Increment time step and go to the next time step right now
  incr timeStep
  return
   
 }


 

 #***** Increment time step and go to the next time step right now
 incr timeStep
 return

} ;# namespace
} ;# proc calcforces {}
 
#*****************
# END
#*****************
