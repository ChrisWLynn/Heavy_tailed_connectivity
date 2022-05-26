# Heavy_tailed_connectivity
Code and data used to perform the analyses in the paper "Heavy-tailed neuronal connectivity arises from Hebbian self-organization" by Christopher W. Lynn, Caroline M. Holmes, and Stephanie E. Palmer.


The neuronal connectome data is included in the following files:

(i) Drosophila_central brain.csv.zip (note that this file must be unzipped): the Drosophila central brain published in Scheffer, L. K. et al. "A connectome and analysis of the adult Drosophila central brain", Elife (2020);

(ii) Drosophila_optic_medulla.csv: the Drosophila optic medulla published in Takemura, S. et al. "A visual motion detection circuit suggested by Drosophila connectomics", Nature (2013);

(iii) Celegans.csv: the C. elegans connectome published in Varshney, L. R., Chen, B. L., Paniagua, E., Hall, D. H. & Chklovskii, D. B. "Structural properties of the Caenorhabditis elegans neuronal network", PLoS Comput. Biol. (2011);

(iv) Platynereis_sensory_motor.csv: the Platynereis sensory-motor circuit published in Randel, N. et al. "Neuronal connectome of a sensory-motor circuit for visual navigation", Elife (2014); and

(v) Mouse_retina.csv: the mouse retina published in Helmstaedter, M. et al. "Connectomic reconstruction of the inner plexiform layer in the mouse retina", Nature (2013).

For each connectome, each row corresponds to an individual synapse (or contact for the mouse retina), with the first, second, and third columns representing the presynaptic neuron, postsynaptic neuron, and synaptic weight, respectively. Note that the synaptic weights are 1 for all connectomes other than the mouse retina (where the synaptic weights have units of squared microns).


The code
