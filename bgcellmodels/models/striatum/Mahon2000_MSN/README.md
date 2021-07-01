# Model Description

This is the Striatal Medium Spiny Neuron (MSN) model presented in Corbit, Whalen
et al (2016) ([ModelDB: 227577](https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=227577)) which is a modification of the model by Mahon, Deniau et al (2000) available at [ModelDB: 150621](https://senselab.med.yale.edu/ModelDB/showmodel.cshtml?model=150621)

# Modifications

The NEURON mod files were downloaded from ModelDB (Mahon) model. Parameters
were adapted according to the modifications decrived in Corbit, Whalen et al. (2016),
i.e. reducing the reversal potential for the leak current from -75 to -90 mV
parameter `seg.el_Leakm` for somatic segments). Mod files were further adapted 
by adding function tables for the channel opening/closing rates.