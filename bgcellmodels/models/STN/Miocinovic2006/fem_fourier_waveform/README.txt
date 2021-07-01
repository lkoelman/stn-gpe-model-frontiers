FEMlab fourier-derived stimulus waveforms for use in NEURON simulations
(normalized to 1V peak-to-peak!!; dt=0.01ms (must be set to that in NEURON))

These files are created with fourier_waveform.m

ffw_soletra_bipo1_130Hz_210us.txt : for bipolar 1 (anode contact 0; cathose contact 2); 130Hz, 210us pulse width
				   based on actual recording from Soletra IPG		
ffw_reconstr_bipo1_encap500_136Hz_90us.txt : for bipolar 1 (anode contact 0; cathose contact 2); 136Hz, 90us pulse width,
					    500um thick encapsulation tissue
				              based on reconstructed waveforms from Intrell II IPG


NOTE: bipo1 and bipo2 files (everything else being the same) are essentially the same	
also mono0,mono1,mono2 and mono3 are the same 

BUT bipos are different from monos --> somethimes they differ by upto 18%