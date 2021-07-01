README for the model files associated with the papers:

1. Khaliq ZM, Gouwens NW, Raman IM (2003) 
The contribution of resurgent sodium current to high-frequency firing in 
Purkinje neurons: an experimental and modeling study. 
J Neurosci 23:4899-912 [PubMed]
2. Raman IM, Bean BP (2001) 
Inactivation and recovery of sodium currents in cerebellar Purkinje 
neurons: evidence for two mechanisms. 
Biophys J 80:729-37

These models were written in NEURON.  Sample use: Download and expand
the archive (zip) file. cd to the directory and run mknrndll (nrnivmodl
on unix).  Then start up the simulation with
nrngui mosinit.hoc

You should see the recreation of graphs that are similar to the model curves
in figure 6 of paper 1 above.  If you want to see a more exact recreation
change the time step to 0.025 (a larger time step of in some cases
1 ms is used in the default.ses file so that it does not take a long time
to reproduce the families of voltage clamp traces).
