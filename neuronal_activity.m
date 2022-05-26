function x  = neuronal_activity(A, b, beta, x0, alpha)
% Inputs: connectivity matrix A, bias vector b, interaction strength beta,
% and vector x0 giving the initial guess at a solution to the self-
% consistent equations. We also include step size of iterations alpha.
%
% Outputs: Solve the self-consistent equation x = tanh(beta(A*x + b)) by
% iterating beginning at x0.

% Convergence threshold:
threshold = 10^(-6);

% Maximum number of steps:
max_steps = 10^6;

% Initialize activities:
x = x0;

% Loop until threshold is reached:
diff = 1;
count = 1;

while diff > threshold
    
    % Change in activities:
    dx = tanh(beta*(A*x + b)) - x;
    diff = max(abs(dx));
    
    % Compute new magnetizations:
    x = x + alpha*dx;
    
    % If threshold has been reached then return:
    if count > max_steps
        error('Maximum number of steps has been reached!');
    end
    
    count = count + 1;
    
end

