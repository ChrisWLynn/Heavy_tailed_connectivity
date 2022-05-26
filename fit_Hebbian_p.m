function [p] = fit_Hebbian_p(Ss, Ps, s_avg, p0)
% Inputs: Observed connection strengths Ss, probabilities of different
% strengths Ps, and average connection strength s_avg. We also take an
% initial guess at the Hebbian probability p0. We note that s_avg is the
% average strength over all possible connections, and therefore may not
% equal sum(Ss.*Ps) if, for example, Ps only runs over positive
% connections.
%
% Outputs: Fit Hebbian probability p in activity-indepent model using the
% analytic connection strength distribution.

% Step size (adjust for convergence):
step_size = 10^(-3);

% Derivative resolution:
dp = 10^(-10);

% Stopping criteria:
dE_dp_min = 10^(-6);
max_steps = 10^4;

% Make sure we only have unique strength values:
Ss_unique = unique(Ss);
Ps_unique = zeros(size(Ss_unique));
num_s = length(Ss_unique);

for i = 1:num_s
    Ps_unique(i) = sum(Ps(Ss == Ss_unique(i)));
end

% Make sure probabilities are normalized to one:
Ps_unique = Ps_unique/sum(Ps_unique);

% Initialize Hebbian probability:
p = p0;

% Loop until max number of steps:
for i = 1:max_steps
    
    % Compute model probabilities:
    Ps_pred = gamma(Ss_unique + s_avg*(1/p - 1))./gamma(Ss_unique + s_avg*(1/p - 1) + 1 + 1/p);
    inds_bad = [find(isnan(Ps_pred)); find(Ps_pred < realmin)];
    if ~isempty(inds_bad)
        c = Ps_pred(min(inds_bad)-1)/Ss_unique(min(inds_bad)-1)^(-(1+1/p));
        Ps_pred(inds_bad) = c*Ss_unique(inds_bad).^(-(1+1/p));
    end
    Ps_pred = Ps_pred/sum(Ps_pred);
    
    % Compute model probabilities slightly away from p (for derivative):
    p_der = p + dp;
    Ps_der = gamma(Ss_unique + s_avg*(1/p_der - 1))./gamma(Ss_unique + s_avg*(1/p_der - 1) + 1 + 1/p_der);
    inds_bad = [find(isnan(Ps_der)); find(Ps_der == 0)];
    if ~isempty(inds_bad)
        c = Ps_der(min(inds_bad)-1)/Ss_unique(min(inds_bad)-1)^(-(1+1/p_der));
        Ps_der(inds_bad) = c*Ss_unique(inds_bad).^(-(1+1/p_der));
    end
    Ps_der = Ps_der/sum(Ps_der);
    
    % Compute error:
    err = sqrt(sum((log(Ps_unique) - log(Ps_pred)).^2));
    err_der = sqrt(sum((log(Ps_unique) - log(Ps_der)).^2));
    
    % Compute gradient of error with respect to Hebbian probability:
    dE_dp = (err_der - err)/dp;
    
    % Stop if derivative is small enough:
    if abs(dE_dp) < dE_dp_min
        return;
    end
    
    % Update Hebbian probability:
%     p = p - step_size*dE_dp;
    p = p - step_size*dE_dp/sqrt(i);
    
    % Stop if p becomes ill-defined:
    if p < 0
        
        p = 0;
        return;
        
    elseif p > 1
        
        p = 1;
        return;
        
    end

end



