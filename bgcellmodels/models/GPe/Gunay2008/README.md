# Porting code from GENESIS

For documentation of syntax and functions used in the original model code,
see [GENESIS manual](http://genesis-sim.org/GENESIS/Hyperdoc/Manual.html).

For format description of GENESIS cell parameter (.p) file, see [readcell documentation](http://www.genesis-sim.org/GENESIS/Hyperdoc/Manual-25.html#readcell).

The celldescription.p files were converted to SWC morphology format in Linux
using [NLMorphologyConverter](http://www.neuronland.org/NL.html) as follows:

```shell
wine start installer.msi # install the program
cd '~/.wine/drive_c/Program Files (x86)/Neuronland/NLMorphologyConverter'
wine NLMorphologyConverter.exe --help
# you can pass file paths in UNIX format to windows executables run via 'wine' command: they will be converted automatically
wine NLMorphologyConverter.exe /home/me/mymorph.p /home/me/mymorph.swc SWC
```

## Units

- GENESIS uses SI units
  - Cm  : `F/m^2`
  - Rm  : `Ohm*m^2`
  - Ra  : `Ohm*m`
  - gbar: `S/m^2`
  - E   : `V`

- NEURON uses following units:
  - Cm : `uF/cm^2  == 1e-6/1e-4 * F/m^2 == 1e-2 * F/m^2`    => x 1e2
  - Rm : `1/gbar `                                          => x 1e4
  - Ra : `Ohm*cm   == 1e-2 * Ohm*m`                         => x 1e2
  - gbar : `S/cm^2 == 1/1e-4 * S/m^2 == 1e4 * S/m^2`        => x 1e-4
  - E : `mV == V*1e-3`                                      => x 1e3

## TODO porting

- DONE: write mechanisms.json
  - [X] find out how Ephys makes SectionLists
      => based on identified SWC types

--------------------------------------------------------------------------------
# GPe cell models

### Gunay (2008)

For model at [ModelDB:114639](https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=114639): 

Gunay, Edgerton, and Jaeger (2008). _Channel Density Distributions Explain Spiking Variability in the Globus Pallidus: A Combined Physiology and Computer Simulation Database Approach._

- main script in `/runs/runsample/setup.g` loads following scripts:

  + `/common/CIPfuncs.g`
      + `/common/simdefaults.g`         -> sets param variables
  
  + `/runs/runsample/readGPparams.g`
      + `/common/GP[i]_default.g`       -> sets param variables
      + `/common/actpars.g`             -> sets param variables
  
  + `/common/make_GP_libary.g`          -> uses param variables
      + `/common/GPchans.g`             -> defines mechanisms using parameters (similiar to .mod file in NEURON)
      + `/common/GPcomps.g`             -> assign biophysical properties and insert channels in compartments defined in .p file


- [X] TODO: look in online [supplemental materials](http://www.jneurosci.org/highwire/filestream/599434/field_highwire_adjunct_files/0/cengiz_p1_manuscript_R2_supptext_djv3_cg2.pdf) to see which exact parameter combinations within ranges in paper Table 1 were used.

    + set `G_<channel>` and `G_mult_<channel>` accordingly
    
    + information in Supplementary Materials:

        + their parameter fitting approach yielded the baseline parameters that
          are specified in their model code (except NaP which is 1.45 in code VS 1.0
          in paper)

        + the multipliers are just used for the different reggions (soma, dend, axon)

        + the ranges given in Table 1 are used for specific experiments


- [X] check if we should use params from Axon hillock or Axon initial segment
      for the axon stub in axonless neuron.

    + DONE: used parameters for axon stub from Hendrickson (2011)
    
    + they are different in the values of active channels
    
    + we are using morphology with single axon equivalent compartment (stub)
      (axonless morphology from Hendrickson (2011))

        + for the axon stub, he uses different parameters than the model
          with axon (see `GP_ax` section in `GP_comps.g`)

        + these parameters are labeled `G_XXX_axon`
          => copy those parameters to Gunay (2008) model

    + if we use morphology `GP1.p` from Gunay (2008):
        
        + look at `GP1.p` file and check distance from soma + diameter for `GP_axHill`
          and `GP_axIS` compartments:

        + => it looks like only the first compartment (first 5 um) is the hillock, 
             and following 10 compartment (next 25 um) are IS.
        
        + => since our stub is 40 um long, this means we need to split it
             and assign both parameter sets.


+ [x] set gradient in `G_Ca_HVA`

    + see end of file `GPcomps.g` : different named sections `GP_dendrite_dXXX_diamYYY`
      get different density
    
    + see files `GP[i].p`: explanation of compartment grouping and associated
         Ca_HVA values
    
    + the mechanisms/parameter combinations defined in `GPcomps.g` are associated
      with specific compartments in this .p file
    
    + the gradients/density mapping to diameter & soma distance is based on the paper
      [Hanson & Smith 2002](https://doi.org/10.1002/cne.10075)

    + easiest way to do this would be to define a new `Ephys.Location` object
      that instantiates itself based on distance from soma and diameter.
      - see function `common/treeutils/path_segments()` for using paths & `distance()`
      - for sec in seclist:
          - get sections on path to soma
          - filter the ones that are in same seclist
          - if diam & distance constraint satisfied + parent/child constraints
            satisfied: return section


### Hendrickson (2011)

For model at [ModelDB:127728](https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=127728)
"Comparison of full and reduced globus pallidus models (Hendrickson 2010)" see:
  
+ main script in `/articleCode/scripts/genesisScripts/GP1axonless_full_synaptic.g`
  First it loads variables from following scripts:
    + `/articleCode/commonGPFull/GP1_defaults.g`
    + `/articleCode/commonGPFull/simdefaults.g`
    + `/articleCode/commonGPFull/actpars.g`

+ The variables are then used in following scripts, in braces: {varname}
  (`GP1axonless_full*.g` -> `make_GP_libary.g` -> ...)
    + `/articleCode/commonGPFull/GP1_axonless.p`
    + `/articleCode/commonGPFull/GPchans.g`
    + `/articleCode/commonGPFull/GPcomps.g`

### Edgerton (2010)

For model at [ModelDB:136315](https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=136315)
"Globus pallidus neuron models with differing dendritic Na channel expression 
(Edgerton et al., 2010)", see: 

+ main script in `/run_example/run_vivo_example.g`, loads scripts:
    + ../common/GP1_constants.g
    + ../common/biophysics/GP1_passive.g
    + ../common/biophysics/GP1_active.g
    + .. see lines with 'getarg' statement for modification of scaling factors/gradients

### Schulthess (2011)

Model at https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=137846

- see main script where following files are loaded:
  + ./paspars.g for passive parameters
  + ../actpars.g for active parameters