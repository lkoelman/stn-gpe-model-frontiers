# Datasets

Datasets for use in basal ganglia simulations.

## File `Brown_Beta_Modulation_Values.txt`

Explanation from PhD thesis Eleanor Dunn:

In addition to the parameters outlined in Chapter3, the model in this study included a modulation of the cortical firing patterns to account for the variation in beta activity observed in
clinical data (Little and Brown, 2012). The variation in beta activity was estimated from experimental
LFP data that was recorded from a single patient with bilateral STN DBS at the
Nuffield Department of Clinical Neurosciences in the University of Oxford. The time-ranging
magnitude of beta-band activity in the experimental LFP signal with DBS off was calculated
by band-pass filtering the LFP data between 10 and 35 Hz, full wave rectifying, and then lowpass
filtering the rectified data with a cutoff frequency of 5 Hz. This provided an estimate of
temporal modulation of the beta frequency band activity of a typical LFP signal in the beta frequency
range with DBS turned off. The LFP modulation was centered at zero by subtracting the
mean, downsampled to 10 Hz, and then added to the direct current input to the cortical soma,
updating every 100 ms during the simulation. This resulted in a modulation of cortical firing
rates and subsequently a variation in the beta-band activity of the LFP over time.