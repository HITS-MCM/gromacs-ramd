namespace eval ::RAMD {

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

    # Assign default values for the parameters not specified in the configuration file
    foreach option [array names defaults] {
        if {! [info exists $option]} {
            set $option $defaults($option)
            puts $ramdfileid [format "RAMD: %25s = %s" $option [expr $$option]]
        } elseif { [info exists $option] } {
            puts $ramdfileid [format "RAMD: %25s = %s" $option [expr $$option]]
        }
    }

    # Check if mandatory parameters are specified in the configuration file
    foreach var $mandatory {
        if {! [info exists $var]} {
            error "RAMD: Mandatory parameter '$var' is not set -- cannot start RAMD"
        } else {
            puts $ramdfileid [format "RAMD: %25s = %s" $var [expr $$var]]
        }
    }

    # Check if 'forceOutFreq' is equal to 1; exit with error if that's the case
    if { $forceOutFreq == 1 } { error "RAMD: ERROR: 'forceOutFreq' parameter may not be 1" } 
   
    # Performed pure RAMD if 'mdSteps' = 0
    if { $mdSteps == 0 } { 
     
        # Check that the number of ramd steps is a multiple of 'forceOutFreq'; exit with error if not
        set r [expr "$ramdSteps % $forceOutFreq"]
        if { $r != 0 } { error "RAMD: ERROR: The number of RAMD steps is not a multiple of 'forceOutFreq'" } 
       
        puts $ramdfileid "RAMD: Pure RAMD simulation is performed" 
        
        # If 'mdSteps' is 0 and "mdStart" is yes, give a warning
        if { $mdStart == "yes" } { 
            puts $ramdfileid "RAMD: WARNING: 'mdStart' has no meaning for pure RAMD simulation; it will be ignored" 
        }
        
        # If 'mdSteps' is 0 and "rMinMd" is set, give a warning
        if { [info exists rMinMd] } {
            puts $ramdfileid "RAMD: WARNING: 'rMinMd' specified while 'mdSteps' is 0"
            puts $ramdfileid "RAMD: WARNING: For combined RAMD-MD simulation 'mdSteps' must be greater than 0"
            puts $ramdfileid "RAMD: WARNING: Ignore 'rMinMd' and perform pure RAMD simulation"
        }
    }
   
    # Perform combined RAMD with MD simulation if 'mdSteps' is not 0 and 'rMinMd' is specified
    if { $mdSteps != 0 } { 
        error "To run ramd setting mdSteps is not supported"
    }
        
    # Make a list of all the atoms on which the force will be applied
    set ramdAtoms {}
    for { set i $firstRamdAtom } { $i <= $lastRamdAtom } { incr i } { lappend ramdAtoms $i }
    puts $ramdfileid "RAMD: Atoms subject to the random acceleration are: $ramdAtoms"
    foreach ramdAtom $ramdAtoms { addatom $ramdAtom }

    # Define a group of the ligand atoms; the force will be applied on the center of mass of this group
    set ramdGroup [ addgroup $ramdAtoms ]
   
    # Define a group containing all protein atoms
    set protAtoms {}
    for { set i $firstProtAtom } { $i <= $lastProtAtom } { incr i } { lappend protAtoms $i }
    foreach protAtom $protAtoms { addatom $protAtom }
    set protGroup [ addgroup $protAtoms ]
   
    # Some variable initialization 
    set timeStep 0; set exitFlag 0; 
    set prevLigCOM "0.0 0.0 0.0"; set prevProtCOM "0.0 0.0 0.0"; set prevDist 0;
   
    # Initialization of simulation flags
    if { $mdSteps == 0 } {
        set ramdFlag 1; set mdFlag 0; set ramdStep 0; set mdStep 0;
    }
    
}; # namespace
 
proc genRandAccel { timeStep } {
namespace eval ::RAMD {

    set pi 3.141592653589793238462643383
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
 
}; # namespace
}; # proc genRandAccel


proc evalWalkDist { timeStep prevLigCOM prevProtCOM currLigCOM currProtCOM } {
namespace eval ::RAMD {
 
    # Compute the relative position of the ligand com with regard to the protein com
    set prevRelLigCOM [ vecsub $prevLigCOM $prevProtCOM ]
    set currRelLigCOM [ vecsub $currLigCOM $currProtCOM ]
    
    # Compute the distance travelled by the ligand com during a ramd or md stint
    set vecWalkDist [vecsub $currRelLigCOM $prevRelLigCOM]
    set walkDist [veclength $vecWalkDist]
   
    return  

}; # namespace
}; # proc evalWalkDist
  
proc calcforces {} {
namespace eval ::RAMD {

    # Terminate NAMD if the ligand has exited from the protein
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
        exit 0
    } 
   
    if { [ array exists coords ] } { array unset coords }
    if { [ array exists totForces ] } { array unset totForces }
   
    # Load coordinates for all the atoms and groups defined
    loadcoords coords

    # Load masses for all the atoms and groups defined
    if {[array exists masses] == 0} { loadmasses masses }

    # Calculate the mass of the ligand
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
   
    # Load total forces from previous time step for all the atoms and groups defined
    loadtotalforces totForces 
   
    # Calculate the position of protein and ligand COM
    set protCOM "$coords($protGroup)"
    set ligCOM "$coords($ramdGroup)"
     
    # Initialize ramd simulation or combined ramd-md simulation that begins with ramd
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
      
        # Initial com positions
        set currLigCOM $ligCOM; set currProtCOM $protCOM  
        puts $ramdfileid "RAMD: $timeStep     >>> LIGAND COM IS: $currLigCOM"
        puts $ramdfileid "RAMD: $timeStep     >>> PROTEIN COM IS: $currProtCOM"
      
        # Evaluate initial distance between ligand com and protein com
        set currDist [veclength [vecsub $currLigCOM $currProtCOM]]
        puts $ramdfileid "RAMD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: DIST = $currDist"
      
        # Generate an initial orientation for the acceleration to be applied when ramd is switched on
        genRandAccel $timeStep
        puts $ramdfileid "RAMD: $timeStep     >>> INITIAL RANDOM DIRECTION: $r :: ||r|| = $lenr"
      
        puts $ramdfileid "RAMD: $timeStep  ***** START WITH $ramdSteps STEPS OF RAMD SIMULATION *****"
       
        # Reset the positions of the ligand and protein COMs and the distance ligand com - protein com
        set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"
      
    # Perform ramd simulation 
    } else { 
    
        # Count ramd steps
        incr ramdStep
        
        # Set flag for writing force output
        if { $forceOutFreq != 0 } { set outputFlag [expr "$ramdStep % $forceOutFreq"] }
        
        # Write force output every 'forceOutFreq' steps 
        if { $outputFlag == 0 } {
            set f [veclength $vecAccel]
            puts $ramdfileid "RAMD FORCE: $timeStep  > LIGAND COM is: $ligCOM\nRAMD FORCE: $timeStep  > PROTEIN COM IS $protCOM\nEXTERNAL FORCE VECTOR (F): $vecAccel; ||F|| = $f\nRAMD FORCE: $timeStep  > EXTERNAL FORCE DIRECTION (r): $r; ||r|| = $lenr"
        }

        # Calculate external and total forces acting on the ligand
        set totLigForce "0 0 0"
        foreach ramdAtom $ramdAtoms {
            set atomForce [ vecscale "$relAtomMass($ramdAtom)" $vecAccel ]
            addforce $ramdAtom $atomForce
            set totAtomForce "$totForces($ramdAtom)"
            set totLigForce [vecadd $totLigForce $totAtomForce]
        }
        
        set totLigForceV [veclength $totLigForce]

        # Write total forces acting on the ligand com
        if { $outputFlag == 0 } {
            puts $ramdfileid "RAMD FORCE: $timeStep  > TOTAL FORCE ON THE LIGAND COM IS: $totLigForce ($totLigForceV)"
        }   

        unset totLigForce;
        unset totLigForceV;
        
        # Set flag for evaluating ramd simulation
        set evalRamdFlag [expr "$ramdStep % $ramdSteps"]
        
        # Evaluate ramd stint
        if { $evalRamdFlag == 0 } {
        
            puts $ramdfileid "RAMD: $timeStep  ***** EVALUATE $ramdSteps RAMD STEPS AT TIMESTEP $timeStep *****"
            
            # com positions
            set currLigCOM $ligCOM; set currProtCOM $protCOM
            
            # Evaluate distance between ligand com and protein com
            set currDist [veclength [vecsub $currLigCOM $currProtCOM]]   
            # Evaluate the change in the distance between the protein and the ligand com during the ramd stint 
            set diffDist [expr "${currDist}-${prevDist}"]
            puts $ramdfileid "RAMD: $timeStep     >>> DISTANCE LIGAND COM - PPROTEIN COM IS: $currDist; IT CHANGED BY $diffDist"
            
            # Check if the ligand has exited the protein
            if { $currDist >= $maxDist } { set exitFlag 1; return }
            
            # Compute the distance travelled by the ligand com during the ramd stint
            evalWalkDist $timeStep $prevLigCOM $prevProtCOM $currLigCOM $currProtCOM
            
            # Evaluate whether a new force direction will be generated
            if { $walkDist <= $rMinRamd } {    
            
                genRandAccel $timeStep
               
                puts $ramdfileid "RAMD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (< $rMinRamd)\nRAMD: $timeStep     >>> CONTINUE WITH $ramdSteps STEPS OF RAMD SIMULATION\nRAMD: $timeStep     >>> CHANGE ACCELERATION DIRECTION TO: $r; ||r|| = $lenr"
                
                # Reset the ramd step count
                set ramdStep 0;

                # Reset the positions of the ligand and protein COMs
                set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"

            } else {

                puts $ramdfileid "RAMD: $timeStep     >>> THE DISTANCE TRAVELLED BY THE LIGAND IS: $walkDist (> $rMinRamd)\nRAMD: $timeStep     >>> CONTINUE WITH $ramdSteps STEPS OF RAMD SIMULATION\nRAMD: $timeStep     >>> KEEP PREVIOUS ACCELERATION DIRECTION: $r; ||r|| = $lenr"
               
                # Reset the positions of the ligand and protein COMs
                set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"
            }  
            # Ensure that the positions of the ligand and protein COMs are reset after the evaluation
            set prevLigCOM "$currLigCOM"; set prevProtCOM "$currProtCOM"; set prevDist "$currDist"
        }
    }
   
    # Increment time step and go to the next time step right now
    incr timeStep
    return

}; # namespace
}; # proc calcforces
 
