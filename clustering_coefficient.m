function cc = clustering_coefficient(A)
% Inputs: Symmetric, binary adjacency matrix A
%
% Output: Clustering coefficient cc, the ratio of the number of triangles
% to the number of triplets in the network

% Check that netwoek is symmetric:
if ~isequal(A, A')
    error('Network is not symmetric!');
end

% Check that network is binary:
if ~all(ismember(unique(A(:)), [0 1]))
    error('Network is not binary!');
end

% Size of network:
N = size(A,1);

% Degrees:
ks = sum(A,2);

% Number of triangles for each node:
triangles = diag(A*triu(A)*A);

% Number of triplets for each node:
triplets = ks(ks > 1).*(ks(ks > 1) - 1)/2;

% Clustering coefficient:
cc = sum(triangles)/sum(triplets);