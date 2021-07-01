# TODO PLAN #

1. first compare balance of currents Otsuka/Gillies/RubinTerman/KhaliqRaman papers and see which currents are responsible for the different features of STN cell dynamics (e.g. bursts, spontaneous pacemaking, ...)
	
	- (plot equilibrium activations variables for each current)
	
	- implement experiments for different firing/physiology features from all papers
		- go through paper Gillies and add experiments to Otsuka test script
		- also add experiments for INaR from DoBean STN/Purkinje & Khaliq Purkinje papers
	
	- run experiment for each feature in each model


2. Construct equivalent reduced model
	
	- reduction method: see principles of computational modeling, Ch. 8 & Ch. 4.3.3
		1. Simplify tree morphology (see Ch. 4.3.3) & p. 85
			- iteratively replace cylinders using technique in Box 3.2 and 3/2 assumption
		2. Tune gbar ration & and geom+electrical properties of equivalent compartment to reproduce important behaviours/experiments of full model
			- tune coupling of two compartments 
				- see cable equation (NEURON book Ch3 p. 26 eq. 3.47 (pdf p. 114), or) Sterrat book eq. 2.24
				- d and Ra determine input resistance for axial current
	
	- pick the balance of currents & implementations that agree most with other STN physiology papers
		- possibly enhance rebound bursts according to Otsuka model to make sure they are robust


3. Add bursting mechanism: removal of slow inactivation of Na currents by IPSPs
	
	- add slowly inactivating Na channels with de-inactivation by IPSPs
		- see refs `DoBean2003/Baufreton2005/Wilson2015`
	
	- then replace the Na current with the more detailed Na channel model of Do/Bean/Khaliq/Raman model
		- but retain the relative magnitude of gnabar and see that you can still reproduce main features of Otsuka & Gillies models
			- potentially include two na mechanisms (one for persistent/transient, one for resurgent)
			- sum of new Na currents must be equal to sum of existing Na currents
	
	- (design experiment to test if resurgent Na currents promote additional patterning unexplained by models without this current)


4. Add synapses
	
	- based on studies synapse distribution
	
	- use synapse mapping procedure Marasco & adapt topology of reduced model to preserve transformations and I/O characteristics of dendritic tree


5. Reduce GPe model
	- see ref

--------------------------------------------------------------------------------
# TODO NEXT #

- Test model reduction code
	
	- [x] Check scaling factors actually used VS paper in RedPurk.hoc
		- same as mentioned in paper: original/equivalent surface (if synapses in cluster: only count sections containing synapse toward original surface)
	
	- [x] See what happens to theoretical Rin in toy model
		- Calculate marasco reduction by hand for simple tree
		- calculated by hand and adapted expressions to make units match
		- in toy tree (secs P/A/B): Rin is conserved witn 0.2% accuracy for toy tree
			- NOTE that cluster root sections are not merged here with Marasco <eq> expressions
	
	- [x] See what happens to theoretical Rin in Gillies & Willshaw model
		- Calculate and compare input resistance of trees in full/reduced model (using algorithm Sterrat Ch. 4)
		- => if trees averaged/RaMERGINGMETHOD=1: they are not equal (see compare_models())
		- => if trees not averaged/RaMERGINGMETHOD=0: they are practically equal (within 1%)
	
	- [x] run other tests and see if all fail
		- => they fail, likely due to different input resistance
	
	- [x] check if mergingYmethod correctly implemented 
		- => NO: i used newri2 for `ri_seq -> diam_seq` (see merge_sec_sequential())
		- however Marasco used newri2 only for `ri_seq` but not for `diam_seq` (see mergingYMethod())
	
	- [ ] Correct spontaneous firing rate through scaling/fitting
		- [x] read Chapt Sterrat parameter fitting
		- [x] compare values post-/pre-reduction
		- [x] test effect of extra dendrites => increases firing rate (due to more I_NaL provided?)
		- => spontaneous firing rate can be adapted by changing cm/gpas/gnaL. However by tuning these parameters you cannot get the firing as low as in the full model (9 Hz)
	
- Implement slow inactivation
	- [ ] get parameters state model slow inactivaton & recovery from papers (see notes below)
		
- See if other currents need to be adjusted to accomodate new current INa_rsg
	- e.g. reduce Ih/HCN or change its parameters
	

--------------------------------------------------------------------------------
# Diary Entries relating to Gillies & Willshaw model

## Diary 03/05/2017 #

- Test incremental reduction methods

	- `FINDING`(`npass=1`, `path_L/path_ri, left_neighbor`): 
		- spontaneous firing faster, but bursts good

	- `FINDING` (`interp_prop=path_L; method=linear_neighbors; npass=1`)
		- burst almost completely identical (Vm both in soma and distally)
		- spontaneous firing only a little faster (`T=84.7`)
	
	- `FINDING` (`npass=7`, trunk sections have no further branching)
		- _plateau protocol_: an accurate burst is maintained until max level of collapsing
		- _rebound protocol_: accurate burst and calcium currents in dendrite
		- _spontaneous protocol_: 
			- no spontaneous firing, likely because baseline INaP during ISI is much lower than in full/spontaneously firing models
			- when I change all `eg.gna_NaL = 8e-6 * 1.3` (i.e. original uniform value * 1.3) from the automatic convex distribution, it works

--------------------------------------------------------------------------------
## Diary 12/04/2017 #

### Observation: identify functional regions based on gbar distribution

- **Method**
	- Run sample.hoc in NEURON GUI
	- use Tools > Model View
	- Heterogeneous parameters > select gbar parameter
	- Check 'SpacePlot' > draw line from start to end point

- **Observations**
	- `gk_KDR`
		- uniform
	- `gk_Kv31`
		- only significant in soma two short sections adjadent to soma
		- peak shape up to soma
	- `gna_Na`
		- only nonzero in soma
	- `gcaT_CaT`
		- zero at soma
		- linearly increasing from next to soma to most distal sections
			- reaches higher value in longest section
	- `gcaL_HVA`
		- stays flat near zero proximally until @ 217 um
			- (halfway short spiny sections, same length along long sections)
		- then it increases linearly to 0.012 in long and 0.0075 in short sections
		- then it stays flat briefly for about 15 um
		- it is slightly higher than near-zero value in soma section
	- `gcaN_HVA`
		- opposite distribution of `gcaL_HCA`:
		- stays flat at zero distally toward proximal 
			- until @ 217 um, halfway short spiny sections, and at same length from soma on long sections
			- i.e. @ 150 um from ends
		- there it starts to increase linearly towards soma
		- in soma itself it is slightly lower than peak
	- `gk_sKCa`
		- nonzero at soma
		- is flat at zero until 240 um from soma (halfway short spiny sections, 2/3 long spiny sections)
		- then there is a step increase to 1e-4 which lasts 70 um
		- in long sections, there is another step increase to twice that value which lasts 36 um

- **Conclusions**
	- 220 um seems to be the points where three of the distributions (`gk_sKCa`, `gcaN_HVA`, `gcaL_HVA`) start to diverge
		- this corresponds to the middle of the long spiny sections

--------------------------------------------------------------------------------
## Diary 20/02/2017 #

### Experiment: run evolutionary optimization Bush & Sejnowski reduction ##

- **Experiment** First optimization
	- => Fittest individual: see CSV files
	- => approx. `soma_diam*0.55, soma_cm*1.45, dend_cm*2, dend_Rm*2, dend_Ra*1.5, dend_diam*0.66, gna_NaL*1.0`
		- _RC factor_ increases everywhere
		- _lambdaAC_ of dendrite increases by factor 3
	- => spiking looks good, although


--------------------------------------------------------------------------------
## Diary 10/02/2017 #

### Experiment: reproduce plateau mechanism in Bush & Sejnowski model ##

- **Experiment** tune bush reduction manually to fix based on plateau insight
	- => interpolation of linear distribution: code verified
	- => INSIGHT: only scaling by area ratios conserves the ratio of gbar in each segment
	- => using the area ratio approach: plateau in dendrite is longer, leading to longer burst

### Experiment: compare relative gbar in original & equivalent models

- **Observation** : compare relative maximum conductances in distal dendritic sections
	- OR = original model
	- EQ = equivalent model
	
	(tree 1, cell 8)
	OR (SThcell[0].dend1[7]) Relative gk_sKCa = 1.0
	OR (SThcell[0].dend1[7]) Relative gcaL_HVA = 79.0
	OR (SThcell[0].dend1[7]) Relative gcaN_HVA = 0.0
	OR (SThcell[0].dend1[7]) Relative gcaT_CaT = 44.0

	(linear_dist, (1,2,4,6,8))
	EQ (spiny) Relative gk_sKCa = 1.0
	EQ (spiny) Relative gcaL_HVA = 62.048173337
	EQ (spiny) Relative gcaN_HVA = 0.0
	EQ (spiny) Relative gcaT_CaT = 46.1754975309

	(linear_dist, (1,3,8))
	EQ (spiny) Relative gk_sKCa = 1.0
	EQ (spiny) Relative gcaL_HVA = 52.5732848281
	EQ (spiny) Relative gcaN_HVA = 0.0
	EQ (spiny) Relative gcaT_CaT = 38.7048076976

	(left_neighbor, (1,3,8))
	EQ (spiny) Relative gk_sKCa = 1.0
	EQ (spiny) Relative gcaL_HVA = 46.166232757
	EQ (spiny) Relative gcaN_HVA = 0.0
	EQ (spiny) Relative gcaT_CaT = 36.6491934546

--------------------------------------------------------------------------------
## Diary 09/02/2017 #

### Experiment: understand rebound/plateau mechanism ##

- **Observation**: comparison conductance distribution in full/reduced model
	
	- _FULL model_
		- gcaT increases slowly and linearly from 0.0012 @ after soma to 0.0042 @ distal dendrite
		- gcaL remains flat @ 0.0001 up to last section (index 8) and there it increases linearly to 0.0075
		- gsKCa remains flat @ 0.0 up to last section (index 8) where it has a step increase @ x=0.25 to 0.0001
	
	- _REDUCED model_
		- gcaT increases slowly and linearly from 0.003 @ after soma to 0.017 @ distal dendrite
		- gcaL in trunk section is flat at 0.0005, smooth section flat at 0.004, in spiny section it increases linearly from 0.0005 to 0.022
		- gsKCa remains flat @ 0.0 up to last sections (spiny) where it has an ~(1-exp(x)) increase to 0.0004

- **Observation**: plot interplay state variables/responsible currents in dendrites
	- When you look at V_dend at the Ca currents + sKCa current:
	- => CaT fast inactivation is slowly disabled (CaT deinactivation) during course of hyperpolarizing pulse
		- this causes I_CaT to spike when hyperpolarizing current is stopped
		- this spike pushes up the activation variable of CaL
		- this triggers the CaL plateau-current
	- => a depolarized plateau is generated in the dendrite (without APs): this propagates to the soma where it causes the burst
	- => CONCLUSION: if plateau doesn't occur in reduced model, it is either not triggered successfully or too much attenuated from dendrites to soma

- **Experiment**: use gbar interpolated via electrotonic length as parameters for linear distributions
	- keep form of proximal/uniform/distal distribution of full model, but scale according to interpolation
	1. implement functions that install a linear distribution (function of electrotonic length)
	2. firt use it to tune bush reduction manually to fix based on plateau insight
	3. adapt optimizer for multiple traces
	4. use optimizer with new functions

--------------------------------------------------------------------------------
## Diary 08/02/2017 #

### Experiment: reproduce firing modes with Bush & Sejnowski reduced model ##

- **Experiment**: disabling area scaling for cm and all gbar
	- => causes much faster spontaneous firing rate
	- => even though in RHS of dV/dt all gbar (numerator) and Cm (denominator) are scaled by the same factor, there is not enough attenuation from dendrites to soma since `lambda ~ sqrt(1/f*C_m*R_a)`
		- i.e. even if the local dV/dt may be the same if you don't scale, to preserve attenuation from dendrites to soma you NEED to scale Cm and therefore also all the gbar

- **Experiment**: reproduce _spontaneous_ activity
	- params: unchanged (gpas & cm scaled by cluster `or_area/eq_area` so tau unchanged)
		- => AHP has same overall shape and repolarization level
		- => f is faster: T=30ms/f=33Hz vs. T=110/f=9
	- params: `RC_factor=2.0/gNaL_factor=0.85` => 
	- params: `RC_factor=1.5/gNaL_factor=0.8` => T=100
	- params: `RC_factor=1.25/gNaL_factor=0.8` => T=86.5
	- params: `RC_factor=1.25/gNaL_factor=0.75` => T=120
	- params: `RC_factor=1.25/gNaL_factor=0.7` => T=inf (no spont. firing)
	- params: `RC_factor=1.15/gNaL_factor=0.75` => T=115
	- params: `RC_factor=1.10/gNaL_factor=0.75` => T=111
	- params: `RC_factor=1.00/gNaL_factor=0.75` => T=105


- **Experiment**: reproduce _plateau_
	- params: `RC_factor=1.5/gNaL_factor=0.8` => T=100
		- => the plateau potential is successfully triggered
		- => the burst has less spikes (4 vs. 6 spikes)
		- => amplitude of these spikes is about the same
	- params: `RC_factor=1.1/gNaL_factor=0.75` => T=111
		- => same result
	- params: `RC_factor=1.0/gNaL_factor=0.75` => T=105
		- => same result

- **Experiment**: reproduce _rebound burst_
	- protocol: depolarize to same Vm level as full model, then release
	- params: `RC_factor=1.0/gNaL_factor=0.75` => T=105
		- => burst contains 4 tightly packed spikes, followed by 4 slower spikes (still much faster than spontaneous rate) that are not really part of the same burst

--------------------------------------------------------------------------------
## Diary 30/01/2017 #


### Experiment: reproduce firing modes with reduced STN model ##

- Experiment: reproduce _spontaneous_ activity
	- params reduced model: gpas @ 75%; cm @ 300% gna_NaL @ 60%
	- => AHP: repolarized to lower level (-62mV vs. -70 mV)
	- => broader AP
		- due to larger RC time constant?
	- => INaT shows second mini-peak because activation peak outlives inactivation peak
		- probably contributes to broader AP
	- => peaks of IKDR and IKv3 are both broader
	- => peaks of Ca currents are broader

- Experiment: reproduce _plateau_
	- params reduced model: gpas @ 75%; cm @ 300% gna_NaL @ 60%
	- => the burst has less spikes (4 vs. 6 spikes)
	- => the spikes have lower amplitude (max 0 mV vs 15 mV) and less re-polarization

- Experiment: reproduce _rebound burst_
	- stimulation current adjusted to -0.15 from -0.25
	- => burst containes less spikes, which are less tightly packed (broader, slower)
	- => the individual current spikes of IKDR and IKv3 are not packed tightly together and riding on a 'plateau' as in full model, but are more widely spaces and return to baseline current level instead

--------------------------------------------------------------------------------
## Diary 29/01/2017 #

- TODO (mail to self)
	- [x] i was only changing g in sec midpoints: need nested. for loop sections/segments
	- [x] try attaching passive dendrite that functions only as leak reservoir
	- [x] Plot on same axis: open, act, inact, with twinx() on left, current with twinx() on right

### Experiment: reduce INaL in all segments of reduced model
- Experiment: reduce `I_NaL` in all segments of reduced model
	- baseline, full model: T = 110 ms (f=9Hz)
	- baseline, reduced model: T = 
	- => gna_NaL @ 90% in all sections => T=21 ms (47 Hz); 
	- => gna_NaL @ 80% in all sections => T=23 ms (43 Hz); INaL,max = -2e-3
	- => gna_NaL @ 70% in all sections => T=27 ms (37 Hz); INaL,max = -1e-3
	- => gna_NaL @ 55% in all sections => T=45 ms (22 Hz); INaL,max = -1.4e-3
	- => gpas @ 200% => T=18 ms
	- => gpas @ 200% => T=20 ms
	- => cm @ 200% => T=27 ms
	- => gpas @ 200%; cm @ 200%; gna_NaL @ 60% => T=103 ms
	- => gpas @ 75%; cm @ 300% gna_NaL @ 60% => T=112 ms
	- => gpas @ 50%; cm @ 300%; gna_NaL @ 60% => T=109 ms
	- => gpas @ 250%; cm @ 250%; gna_NaL @ 60% => T=inf: no spontaneous firing
	- => gpas @ 300%; cm @ 300%; gna_NaL @ 60% => T=inf: no spontaneous firing

- Experiment: attach passive tree
	- => with passive copy of trunk0 tree => T=24 ms; INaL,max = -2.5e-3

--------------------------------------------------------------------------------
## Diary 19/01/2017 #

- Try today
	- [x] take another look as scaling g with surface method (try direct surface scaling approach)
		=> same result
	- [x] only lower gNaL of soma
	- [x] plot (in) activation variables and study role in AHP => control AHP duration by setting balance of g

### Experiment & Insight: role of ionic currents (gbar) ##

- [x] plot (in) activation variables and study role in AHP => control AHP duration by setting balance of g
	
	- increase `g{sKCa}`
		+ sKCa is (re)polarizing
		+ sKCa prolongs AHP and prevent bootstrapping/run-away of CaT
		+ SKCa competes with CaT to suppress CaL bootstrapping/runaway
		+ delays the takeover of the membrane potential by the INaP (persistent sodium current)
	
	- dechrease `Ih`
		+ hyperpolarization-activated depolarizing current
			- Reversal potential is -5 mV in our model -> hyperpolarizes during spike high phase
		+ speeds transition into the activation range of INaP
		+ activates slowly at hyperpolarized level
    	+ together with CaT, causes initial depolarization after sustained hyperpolarization
    
    - decrease `g{CaT}`
    	+ activated at low voltage threshold
    	+ - CaT bootstraps CaL (activates at hyperpolarized V and inactivates at depolarized V)

--------------------------------------------------------------------------------
## Diary 11/01/2017 #

### Experiments ##
Test spontaneous firing
	- [x] simulated full model
		=> rgular spiking @ 10 Hz
	- [x] Simulated using average axial resistance (Marasco method)
		=> regular spiking @ 39 Hz
	- [x] Simulated using conservation of Rin (no averaging of trees)

Test rebound burst simulator protocol
	- [x] simulated full model
	- [x] Simulated using average axial resistance (Marasco method)
	- [x] Simulated using conservation of Rin (no averaging of trees)

Test generation of plateau potential
	- [x] simulated full model
	- [x] Simulated using average axial resistance (Marasco method)
	- [x] Simulated using conservation of Rin (no averaging of trees)

### Fitting/Rescaling ###

- to reduce firing rate:
	- scaling L (surface) does not work
	- scaling Cm and Rm/gpas does not work
	- applying a surface scaling factor like in Marasco paper does not work

- conclusion: fit passive RC circuit model to full model as suggested in Sterrat book