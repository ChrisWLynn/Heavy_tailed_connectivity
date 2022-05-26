# Heavy_tailed_connectivity
Code and data used to perform the analyses in the paper "Heavy-tailed neuronal connectivity arises from Hebbian self-organization" by Christopher W. Lynn, Caroline M. Holmes, and Stephanie E. Palmer.

The neuronal connectomes are included in the following files:

(i) Drosophila_central brain.csv.zip (note that this file must be unzipped): the Drosophila central brain published in Scheffer, L. K. et al. "A connectome and analysis of the adult Drosophila central brain", Elife (2020);

(ii) Drosophila_optic_medulla.csv: the Drosophila optic medulla published in Takemura, S. et al. "A visual motion detection circuit suggested by Drosophila connectomics", Nature (2013);

(iii) Celegans.csv: the C. elegans connectome published in Varshney, L. R., Chen, B. L., Paniagua, E., Hall, D. H. & Chklovskii, D. B. "Structural properties of the Caenorhabditis elegans neuronal network", PLoS Comput. Biol. (2011);

(iv) Platynereis_sensory_motor.csv: the Platynereis sensory-motor circuit published in Randel, N. et al. "Neuronal connectome of a sensory-motor circuit for visual navigation", Elife (2014); and

(v) Mouse_retina.csv: the mouse retina published in Helmstaedter, M. et al. "Connectomic reconstruction of the inner plexiform layer in the mouse retina", Nature (2013).

For each connectome, each row corresponds to an individual synapse (or contact for the mouse retina), with the first, second, and third columns representing the presynaptic neuron, postsynaptic neuron, and synaptic strength, respectively. Note that the synaptic strengths are 1 for all connectomes other than the mouse retina (where the synaptic strengths have units of squared microns).

To analyze the connectomes, use "analyze_connectomes.m". To plot the distributions of the connection strengths, use "plot_connection_strengths.m". To simulate the activity-independent model, use either "activity_independent_directed.m" (for directed networks) or "activity_independent_symmetric.m" (for symmetric networks). To simulate the activity-dependent model, use "activity_dependent.m".

The above scripts use the following helper functions:

"clustering_coefficient.m" computes the clustering coefficient for a given symmetric connectivity matrix.

"fit_Hebbian_p.m" fits the Hebbian probability p (using the analytic connection strength distribution for the activity-independent model) to a given connection strength distribution.

"neuronal_activity.m" computes the average activities of model neurons obeying the self-consistency equation x = tanh[beta*(Ax + b)], which is used in the activity-dependent model.
