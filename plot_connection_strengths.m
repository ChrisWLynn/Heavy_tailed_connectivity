%% Script to plot distributions of connection strengths

% Index of connectome to analyze:
% ind = 1; % Drosophila central brain
% ind = 2; % Drosophila optic medulla
% ind = 3; % C elegans
% ind = 4; % Platynereis sensory-motor circuit
ind = 5; % Mouse retina

% Connectomes:
filenames = {'Drosophila_central_brain', 'Drosophila_optic_medulla', 'Celegans',...
    'Platynereis_sensory_motor', 'Mouse_retina'};

% Load connectome:
synapses = readmatrix([filenames{ind}, '.csv']);

% Number of neurons:
N = max(synapses(:, [1 2]), [], [1 2]);

% Directed connectivity matrix:
A = full(sparse(synapses(:,1), synapses(:,2), synapses(:,3), N, N));

% Connection strengths:
conn_strengths = A(A > 0);

% Compute connection strength distribution:
bins = 0:ceil(max(conn_strengths));
s_unique = zeros(1, length(bins) - 1);
Ps = zeros(1, length(bins) - 1);

for i = 1:(length(bins) - 1)

    inds_bin = logical((conn_strengths > bins(i)).*(conn_strengths <= bins(i+1)));
    
    if sum(inds_bin) > 0
        s_unique(i) = mean(conn_strengths(inds_bin));
        Ps(i) = sum(inds_bin)/length(conn_strengths);
    end
end

s_unique = s_unique(Ps > 0);
Ps = Ps(Ps > 0);

% Fit Hebbian probability for activity-independent model:
p = fit_Hebbian_p(s_unique', Ps, 1, .5);

%% Plot connection strength distribution:

% Plot options:
colors = [0 77 128; 181 23 0; 1 113 0; 242 112 0; 120 0 150; 0 168 157; 203 41 123; 0 0 0]/255;
marker_alpha = .8;
line_width = 2;
axes_width = 1;
font_size = 17;
marker_size = 100;
x_limits = [1, 2*10^3; 1, 5*10^2; 1, 5*10^1; 1, 5*10^1; 10^(-1), 2*10^2];
y_limits = [10^(-7), 1; 10^(-4), 1; 2*10^(-4), 1; 2*10^(-3), 1; 10^(-5), 1];

figure;
hold on;

% Plot connectivity strength distribution:
scatter(s_unique, Ps, marker_size, marker_alpha*colors(ind,:) + 1-marker_alpha, 'filled', 'MarkerFaceAlpha', 1);

% Plot analytic activity-independent distribution:
Ps_Hebb = gamma(s_unique + 1/p - 1)./gamma(s_unique + 2/p);
inds_bad = [find(isnan(Ps_Hebb)), find(Ps_Hebb < realmin)];
if ~isempty(inds_bad)
    c = Ps_Hebb(min(inds_bad)-1)/s_unique(min(inds_bad)-1)^(-(1+1/p));
    Ps_Hebb(inds_bad) = c*s_unique(inds_bad).^(-(1+1/p));
end
Ps_Hebb = Ps_Hebb/sum(Ps_Hebb);
plot(s_unique, Ps_Hebb, '-', 'Color', colors(end,:), 'LineWidth', line_width + 1.5)

% Plot appearance:
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.XLim = x_limits(ind,:);
ax.YLim = y_limits(ind,:);
ax.LineWidth = axes_width;
ax.FontSize = font_size;
ax.TickDir = 'out';

pbaspect([1 1 1])
box('off')
hold off;

