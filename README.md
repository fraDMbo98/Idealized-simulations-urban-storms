# Idealized simulations of urban storms

Code and input files used to perform and analyse the idealized WRF simulations of
deep convection over an urban area presented in:

> De Martin, F., Rozoff, C., Zonato, A., Alessandrini, S., and Di Sabatino, S.:
> The Dual Effect of Urban Areas on Supercell Storms, EGUsphere [preprint],
> <https://doi.org/10.5194/egusphere-2026-815>, 2026.

Archived version: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23110683.svg)](https://doi.org/10.5281/zenodo.23110683)

## Repository structure

```
WRF/
├── namelist.input             WRF namelist
├── input_sounding             initial thermodynamic and wind profile
├── surface_file               soil/vegetation types, initial soil moisture, orography parameters
├── city_file                  city position and extent
├── bubble_file                optional warm bubble (position, amplitude, time)
├── URBPARM_LCZ.TBL            urban parameters for the local climate zones
├── MPTABLE.TBL                Noah-MP parameters (NATURAL = 12)
├── solve_em.F                 delayed warm bubble, passive tracers, time-averaged fields
├── module_initialize_ideal.F  idealized case urban_xy
├── module_sf_bep_bem.F        BEP+BEM urban canopy model (development version)
├── module_sf_noahmpdrv.F      Noah-MP driver (control), standard urban–rural coupling
└── WRF_nothermo/              files used only for the no_thermo experiment
    ├── module_sf_noahmpdrv.F  Noah-MP driver, urban thermal and radiative fluxes replaced by rural ones
    └── bl_ysu.f90             YSU PBL scheme without urban heat and moisture sources
SPOUT/                         revised SPOUT code for updraft identification and tracking
jupyter/                       notebooks used to process the output and produce the figures
```

## How to run the simulations

1. Download WRF **v4.6.0** from <https://github.com/wrf-model/WRF>.
2. Copy the source files into the WRF tree, replacing the original ones:
   - `solve_em.F`, `module_initialize_ideal.F` → `dyn_em/`
   - `module_sf_bep_bem.F` → `phys/`
   - `module_sf_noahmpdrv.F` → `phys/noahmp/drivers/wrf/`
   - only for the *no_thermo* experiment, replace the previous files with those in `WRF_nothermo/`:
     `module_sf_noahmpdrv.F` → `phys/noahmp/drivers/wrf/` and
     `bl_ysu.f90` → `phys/physics_mmm/`
3. Compile the idealized case (`./configure`, then `./compile em_<case>`).
4. Copy `namelist.input`, `input_sounding`, `surface_file`, `city_file`,
   `bubble_file` and the two `.TBL` tables into the run directory, then run
   `ideal.exe` followed by `wrf.exe`.

The control and *no_thermo* experiments use the same configuration and differ
only in the files contained in `WRF_nothermo/`. In the *no_thermo*
experiment, the urban canopy still acts mechanically (building drag and
volume), while radiation, albedo, emissivity and the sensible, latent and
ground heat fluxes of urban grid cells are taken from the rural surface.

## Notes

- `bl_ysu.f90` is the preprocessed file from the build directory: running
  `./clean` restores the original `bl_ysu.F90`, so copy it again before recompiling.
- Model outputs are not included because of their size and are available from
  the corresponding author upon request.

## Contact

Francesco De Martin — francescodemartin1998@gmail.com
