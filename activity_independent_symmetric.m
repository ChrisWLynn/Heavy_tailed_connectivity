% Script to simulate symmetric activity-independent model

% Number of neurons:
N = 1000;

% Number of possible connections:
E = N*(N-1)/2;

% Average connection strength:
s_avg = 1;

% Number of undirected synapses:
num_syn = round(s_avg*E);

% Probability of Hebbian growth:
p = 0.5;

% Number of data samples:
num_samples = 100;

% Number of network updates per data sample:
num_updatesPerSample = 100;

% Number of parallel connection prunes per network update:
num_prunesPerUpdate = ceil(E/num_updatesPerSample);

% Number of network updates for burn-in:
num_updatesBurn = 50*num_updatesPerSample;

% Quantities to compute:
s_values = cell(num_samples, 1);
s_counts = cell(num_samples, 1);
density = zeros(num_samples, 1);
heterogeneity = zeros(num_samples, 1);
clustering = zeros(num_samples, 1);
        
% Initialize network with synapses placed randomly:
edge_inds = find(triu(ones(N), 1));
[I, J] = ind2sub([N, N], edge_inds);
inds_sample = randsample(E, num_syn, true);
A = full(sparse(I(inds_sample), J(inds_sample), ones(num_syn,1), N, N));
A = A + A';

% Loop over burn-in steps:
for i = 1:num_updatesBurn
    
    % Pick connections to prune:
    inds_remove = randsample(E, num_prunesPerUpdate, false);
    
    % Prune connections:
    s_temp = sum(A(edge_inds(inds_remove)));
    A = A - full(sparse([I(inds_remove); J(inds_remove)], [J(inds_remove); I(inds_remove)], [A(edge_inds(inds_remove)); A(edge_inds(inds_remove))], N, N));
    
    % Number of Hebbian updates:
    s_Hebb = binornd(s_temp, p);
    
    % Pick connections to increase with Hebbian growth:
    inds_inc_Hebb = randsample(E, s_Hebb, true, A(edge_inds) + realmin);
    
    % Pick connections to increase randomly:
    inds_inc_rand = randsample(E, s_temp - s_Hebb, true);
    
    % Update connection strengths:
    inds_inc = [inds_inc_Hebb; inds_inc_rand];
    A = A + full(sparse([I(inds_inc); J(inds_inc)], [J(inds_inc); I(inds_inc)], ones(2*length(inds_inc), 1), N, N));
    
end

% Loop over samples:
for i = 1:num_samples
    
    % Loop over network updates per sample:
    for j = 1:num_updatesPerSample
        
        % Pick connections to prune:
        inds_remove = randsample(E, num_prunesPerUpdate, false);
        
        % Prune connections:
        s_temp = sum(A(edge_inds(inds_remove)));
        A = A - full(sparse([I(inds_remove); J(inds_remove)], [J(inds_remove); I(inds_remove)], [A(edge_inds(inds_remove)); A(edge_inds(inds_remove))], N, N));
        
        % Number of Hebbian updates:
        s_Hebb = binornd(s_temp, p);
        
        % Pick connections to increase with Hebbian growth:
        inds_inc_Hebb = randsample(E, s_Hebb, true, A(edge_inds) + realmin);
        
        % Pick connections to increase randomly:
        inds_inc_rand = randsample(E, s_temp - s_Hebb, true);
        
        % Update connection strengths:
        inds_inc = [inds_inc_Hebb; inds_inc_rand];
        A = A + full(sparse([I(inds_inc); J(inds_inc)], [J(inds_inc); I(inds_inc)], ones(2*length(inds_inc), 1), N, N));
        
    end
    
    % Compute quantities:
    
    % Connection strengths:
    conn_strengths = A(edge_inds);
    conn_strengths = conn_strengths(conn_strengths > 0);
    
    % Counts of connections with different weights:
    s_values{i} = unique(conn_strengths)';
    s_counts{i} = histcounts(conn_strengths, [s_values{i}, inf], 'Normalization', 'count');
    
    % Density, clustering coefficient, and heterogeneity:
    density(i) = length(conn_strengths)/E;
    clustering(i) = clustering_coefficient(double(A > 0));
    heterogeneity(i) = 1/2*sum(abs((s_values{i}' - s_values{i}).*...
        ((s_counts{i}'/sum(s_counts{i}))*(s_counts{i}/sum(s_counts{i})))), [1 2])/...
        mean(conn_strengths);
    
end