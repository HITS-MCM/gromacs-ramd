# Extract coord from checkpoint file and trajectory

```
cd /workspaces/gromacs-ramd/build/src/programs/mdrun/tests/Testing/Temporary
/workspaces/gromacs-ramd/build/bin/gmx dump -cp RAMDTest_RAMD_connected_ligands_state.cpt > RAMDTest_RAMD_connected_ligands_state.cpt.dump
/workspaces/gromacs-ramd/build/bin/gmx dump -f RAMDTest_RAMD_connected_ligands.trr > RAMDTest_RAMD_connected_ligands.trr.dump
```
