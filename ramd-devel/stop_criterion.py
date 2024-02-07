# The binding site residues (bind_res) are detected during the conventional MD simulations run before RAMD. 
# The user adds the bind_res of proteins 1 and 2, by creating two new groups in the gromacs index file (e.g., calling them bind_res_1 and bind_res_2).

# 2. Storage of contacts during RAMD trajectories (ramd_contacts).

# The COM-COM distances between bind_res_1 and bind_res_2 are computed during the RAMD trajectory.

# Below is the piece of MDAnalysis code used to compute and store the interactions (i.e. COM-COM distances):

# First, different classes of atoms are defined based on the type of interaction:

at_aromatic="((resname PHE TRP TYR HIS HIE HID HE2) and (name CZ* CD* CE* CG* CH* NE* ND*))"
at_positive="((resname ARG LYS ) and (name NH* NZ)) or ((resname HIP ) and (name HD HE))"   #please note that for resname HIP also, name ND NE can be used
at_negative=" ((resname ASP GLU) and (name OE* OD*))"
at_Hdon="(resname TYR SER LYS GLN ARG HIS ASN THR and (not backbone) and (name O* N*))"
at_Hacc=" (resname GLU ASP GLN and (not backbone) and (name O*))"
at_sulfur="(protein and (name S*))"
at_hydrophob="(protein and (name C* S*) and (not (name CG and resname ASN ASP)) and (not (name CD and resname GLU GLN ARG)) and (not (name CZ and resname TYR ARG)) and (not (name CE and resname LYS)) and (not (name CB and resname SER THR)) and (not backbone))"
at_BB_positive="backbone and name N"
at_BB_negative="backbone and name O"

#define bind_res_1 and bind_res_2
resi_lig = bind_res_1
resi_target = bind_res_2

# Create lists to append the group of atoms for a specific interaction, for both resi_lig and resi_target
  u_CA_lig_IP = []
  u_CA_lig_IN = []
  u_CA_lig_HD = []
  u_CA_lig_HA = []
  u_CA_lig_AR = []
  u_CA_lig_HY = []
  u_CA_lig_BB = []
  u_CA_lig_BB_P = []
  u_CA_lig_BB_N = []
  u_CA_lig_C = []
  u_CA_target_IP = []
  u_CA_target_IN = []
  u_CA_target_HD = []
  u_CA_target_HA = []
  u_CA_target_AR = []
  u_CA_target_HY = []
  u_CA_target_BB = []
  u_CA_target_BB_P = []
  u_CA_target_BB_N = []
  u_CA_target_C = []

  for l in range(resi_lig[0], resi_lig[1]+1 ):
        u_CA_lig_IP.append(u.select_atoms("resid " +str(l)+" and "+at_positive , updating=True))
        u_CA_lig_IN.append(u.select_atoms("resid " +str(l)+" and "+at_negative , updating=True))
        u_CA_lig_HD.append(u.select_atoms("resid " +str(l)+" and "+at_Hdon, updating=True))
        u_CA_lig_HA.append(u.select_atoms("resid " +str(l)+" and "+at_Hacc, updating=True))
        u_CA_lig_AR.append(u.select_atoms("resid " +str(l)+" and "+at_aromatic, updating=True))
        u_CA_lig_HY.append(u.select_atoms("resid " +str(l)+" and "+at_hydrophob, updating=True))
        u_CA_lig_BB.append(u.select_atoms("resid " +str(l)+" and backbone", updating=True))
        u_CA_lig_BB_N.append(u.select_atoms("resid " +str(l)+" and backbone and name O*", updating=True))
        u_CA_lig_BB_P.append(u.select_atoms("resid " +str(l)+" and backbone and name N*", updating=True))
        u_CA_lig_C.append(u.select_atoms("resid " +str(l)+" and name C* and not backbone", updating=True))

  for t in range(resi_target[0], resi_target[1]+1):
        u_CA_target_IP.append(u.select_atoms("resid " +str(t)+" and "+at_positive , updating=True))
        u_CA_target_IN.append(u.select_atoms("resid " +str(t)+" and "+at_negative , updating=True))
        u_CA_target_HD.append(u.select_atoms("resid " +str(t)+" and "+at_Hdon, updating=True))
        u_CA_target_HA.append(u.select_atoms("resid " +str(t)+" and "+at_Hacc, updating=True))
        u_CA_target_AR.append(u.select_atoms("resid " +str(t)+" and "+at_aromatic, updating=True))
        u_CA_target_HY.append(u.select_atoms("resid " +str(t)+" and "+at_hydrophob, updating=True))
        u_CA_target_BB.append(u.select_atoms("resid " +str(t)+" and backbone", updating=True))
        u_CA_target_BB_N.append(u.select_atoms("resid " +str(t)+" and backbone and name O*", updating=True))
        u_CA_target_BB_P.append(u.select_atoms("resid " +str(t)+" and backbone and name N*", updating=True))
        u_CA_target_C.append(u.select_atoms("resid " +str(t)+" and name C* and not backbone", updating=True))

# Calculate the interactions for every pair of residues

for l in range(0, len(u_CA_lig_IP)):
        lt = []
        for t in range(0, len(u_CA_target_IP)):
             # possible polar bond
             if len(u_CA_target_C[t]) == 0 or len(u_CA_lig_C[l]) == 0:
                  d = int(10* LA.norm (u_CA_lig_BB[l].center_of_mass()- u_CA_target_BB[t].center_of_mass()))
             else:
                  d = int(10* LA.norm (u_CA_lig_C[l].center_of_mass()- u_CA_target_C[t].center_of_mass()))
             tt = []
             tt.append(str(d)+"--")
             if d < 150:
               if len(u_CA_lig_IP[l]) > 0 and len(u_CA_target_IN[t]) > 0:
                  d1 = int(10* LA.norm (u_CA_lig_IP[l].center_of_mass()- u_CA_target_IN[t].center_of_mass()))
                  tt.append("IP-IN"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_IN[l]) > 0 and len(u_CA_target_IP[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_IN[l].center_of_mass()- u_CA_target_IP[t].center_of_mass()))
                  tt.append("IN-IP"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HD[l]) > 0 and len(u_CA_target_HA[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HD[l].center_of_mass()- u_CA_target_HA[t].center_of_mass()))
                  tt.append("HD-HA"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HA[l]) > 0 and len(u_CA_target_HD[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HA[l].center_of_mass()- u_CA_target_HD[t].center_of_mass()))
                  tt.append("HA-HD"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_AR[l]) > 0 and len(u_CA_target_AR[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_AR[l].center_of_mass()- u_CA_target_AR[t].center_of_mass()))
                  tt.append("AR-AR"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HY[l]) > 0 and len(u_CA_target_HY[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HY[l].center_of_mass()- u_CA_target_HY[t].center_of_mass()))
                  tt.append("HY-HY"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_BB_P[l]) > 0 and len(u_CA_target_HA[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_BB_P[l].center_of_mass()- u_CA_target_HA[t].center_of_mass()))
                  tt.append("BB-HA"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_BB_N[l]) > 0 and len(u_CA_target_HD[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_BB_N[l].center_of_mass()- u_CA_target_HD[t].center_of_mass()))
                  tt.append("BB-HD"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HD[l]) > 0 and len(u_CA_target_BB_N[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HD[l].center_of_mass()- u_CA_target_BB_N[t].center_of_mass()))
                  tt.append("HD-BB"+str(d)+"-"+str(d1))
                  d = min(d1,d)
               if len(u_CA_lig_HA[l]) > 0 and len(u_CA_target_BB_P[t]) > 0:
                  d1 = int(10* LA.norm(u_CA_lig_HA[l].center_of_mass()- u_CA_target_BB_P[t].center_of_mass()))
                  tt.append("HA-BB"+str(d)+"-"+str(d1))
                  d = min(d1,d)
             lt.append(d)
             
# 3. Residence time computation

# The frame after the time when the mean of ramd_contacts [bind_res] is > 5.5 Angstrom, is defined as the residence time.
# If the condition is not satisfied; terminate the simulation with the message “No dissociation, simulation terminated early”.

# Below is the Python function used, where:

# k= matrix containing the COM-COM distances between bind_res_1 and bind_res_2 (computed at step 2) over time
# threshold=5.5
# tr = the trajectory 

def fun_resi_first(k,threshold,tr):
    k1 = k.time.values[k.values.mean(axis=1) > threshold*10 ]
    if len(k1) > 1:
        time = k1[1]
    else:
        print ("Traj",tr,"The first unbinding event was not found: no dissociation? " )
        time = k.time.values[-1]
    return(time)
