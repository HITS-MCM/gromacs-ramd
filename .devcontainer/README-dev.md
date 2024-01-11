Load development modules on HPC cluster

```bash
ml purge
ml unuse $MODULEPATH
ml use /hits/sw/its/doserbd/cascade/modules/all
ml foss/2022a 
ml CUDA/11.7.0
```
