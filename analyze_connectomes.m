%% Script to compute values of interest for a given connectome

% Connectome to analyze:
% filename = 'Drosophila_central_brain';
% filename = 'Drosophila_optic_medulla';
filename = 'Celegans';
% filename = 'Platynereis_sensory_motor';
% filename = 'Mouse_retina';

% Load connectome:
synapses = readmatrix([filename, '.csv']);

% Number of neurons:
N = max(synapses(:, [1 2]), [], [1 2]);

%% Compute quantities for directed connectome:

% Directed connectivity matrix:
A = full(sparse(synapses(:,1), synapses(:,2), synapses(:,3), N, N));

% Total synaptic strength (or contact area for mouse retina):
S = sum(A(:));

% Average connection strength (or contact area for mouse retina):
s_avg = S/(N*(N-1));

% Number of synapses (or contacts for mouse retina):
num_syn = size(synapses,1);

% Number of connections:
num_conn = sum(A(:) > 0);

% Connection density:
density = num_conn/(N*(N-1));

% Connection heterogeneity:
conn_strengths = A(A > 0);
bins = 0:ceil(max(conn_strengths));
bin_values = zeros(1, length(bins) - 1);
bin_counts = zeros(1, length(bins) - 1);
for i = 1:(length(bins) - 1)
    inds = logical((conn_strengths > bins(i)).*(conn_strengths <= bins(i+1)));
    if sum(inds) > 0
        bin_values(i) = mean(conn_strengths(inds));
        bin_counts(i) = sum(inds);
    end
end
bin_values = bin_values(bin_counts > 0);
bin_counts = bin_counts(bin_counts > 0);
heterogeneity = 1/2*sum(abs(bin_values' - bin_values).*((bin_counts'/sum(bin_counts))*(bin_counts/sum(bin_counts))), [1 2])/...
    sum(bin_values.*bin_counts/sum(bin_counts));

%% Compute quantities for symmetrized connectome:

% Symmetrized connectome:
A_sym = A + A';

% Number of connections:
num_conn_sym = sum(A_sym(:) > 0)/2;

% Connection density:
density_sym = num_conn_sym/(N*(N-1)/2);

% Connection heterogeneity:
conn_strengths_sym = A_sym(triu(A_sym, 1) > 0);
bins = 0:ceil(max(conn_strengths_sym));
bin_values = zeros(1, length(bins) - 1);
bin_counts = zeros(1, length(bins) - 1);
for i = 1:(length(bins) - 1)
    inds = logical((conn_strengths_sym > bins(i)).*(conn_strengths_sym <= bins(i+1)));
    if sum(inds) > 0
        bin_values(i) = mean(conn_strengths_sym(inds));
        bin_counts(i) = sum(inds);
    end
end
bin_values = bin_values(bin_counts > 0);
bin_counts = bin_counts(bin_counts > 0);
heterogeneity_sym = 1/2*sum(abs(bin_values' - bin_values).*((bin_counts'/sum(bin_counts))*(bin_counts/sum(bin_counts))), [1 2])/...
    sum(bin_values.*bin_counts/sum(bin_counts));

% Clustering coefficient:
A_sym_binary = double(A_sym > 0);
clustering_sym = clustering_coefficient(A_sym_binary);

%% Compute quantities for synapse-shuffled network:

% List of possible undirected connections:
[I, J] = ind2sub([N, N], find(triu(ones(N),1)));
num_conn_possible = length(I);

% Generate random synapse-shuffled network:
inds = randsample(num_conn_possible, num_syn, true);
A_randSyn = full(sparse(I(inds), J(inds), synapses(:,3), N, N));
A_randSyn = A_randSyn + A_randSyn';

% Number of connections:
num_conn_randSyn = sum(A_randSyn(:) > 0)/2;

% Connection density:
density_randSyn = num_conn_randSyn/(N*(N-1)/2);

% Connection heterogeneity:
conn_strengths_randSyn = A_randSyn(triu(A_randSyn, 1) > 0);
bins = 0:ceil(max(conn_strengths_randSyn));
bin_values = zeros(1, length(bins) - 1);
bin_counts = zeros(1, length(bins) - 1);
for i = 1:(length(bins) - 1)
    inds = logical((conn_strengths_randSyn > bins(i)).*(conn_strengths_randSyn <= bins(i+1)));
    if sum(inds) > 0
        bin_values(i) = mean(conn_strengths_randSyn(inds));
        bin_counts(i) = sum(inds);
    end
end
bin_values = bin_values(bin_counts > 0);
bin_counts = bin_counts(bin_counts > 0);
heterogeneity_randSyn = 1/2*sum(abs(bin_values' - bin_values).*((bin_counts'/sum(bin_counts))*(bin_counts/sum(bin_counts))), [1 2])/...
    sum(bin_values.*bin_counts/sum(bin_counts));

%% Compute quantities for connection-shuffled network:

% Generate random connection-shuffled network:
inds = randsample(num_conn_possible, num_conn_sym, false);
A_randConn = full(sparse(I(inds), J(inds), conn_strengths_sym, N, N));
A_randConn = A_randConn + A_randConn';

% Clustering coefficient:
A_randConn_binary = double(A_randConn > 0);
clustering_randConn = clustering_coefficient(A_randConn_binary);

