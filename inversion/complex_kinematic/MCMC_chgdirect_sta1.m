clc; clear; close all;

% === Load observed frequency data ===
data = load('../../data/max_freqs_chg1.txt');

Time_obs = data(:,1)';   % Observation time (s)
f_obs = data(:,2);       % Observed frequencies (Hz)
w_obs = ones(size(f_obs));  % Default weights for each observation point

% --- Model parameters ---
% f0 = m_pred(1);             % Source frequency (Hz)
% v_mag = m_pred(2);          % Horizontal rupture velocity (m/s)
% L1 = m_pred(3);             % Perpendicular distance from source to station (m)
% L2 = m_pred(4);             % Perpendicular distance from source to station (m)
% t0 = m_pred(5);             % Source activation time (s)
% theta_end_deg = m_pred(6);  % Direction angle (degrees)
% R = m_pred(7);              % Arc radius (m)

% === Prior distribution (Gaussian) ===
prior_mu  = [25, 100, 500, 500, 20, 135, 200];
prior_std = [0.1, 0.1, 10, 10, 0.1, 0.1, 0.1];

% Log-prior (Gaussian form)
log_prior = @(m) -0.5 * sum(((m - prior_mu)./prior_std).^2);

% === Likelihood function ===
log_likelihood = @(m) ...
    -0.5 * sum(w_obs .* (f_obs - doppler_fwd_chgdirection(m, Time_obs, 340)).^2, 'omitnan');

% Posterior (unnormalized)
log_posterior = @(m) log_prior(m) + log_likelihood(m);

% === MCMC Setup ===
n_iter = 2000;          % Total MCMC iterations
burn_in = n_iter * 0.4;  % Burn-in proportion
m_dim = 7;               % Number of parameters

chain = zeros(n_iter, m_dim);
chain(1,:) = prior_mu;   % Initialize chain with prior mean

% Adaptive step size parameters
step_size = [1, 1, 1, 1, 1, 1, 1] / 100;  
step_size_min = step_size / 100;
step_size_max = step_size * 10;

accept_count = 0;
adapt_interval = 100;    % Step size adjustment interval (not used here explicitly)

% === Main MCMC Loop ===
for i = 2:n_iter
    m_current = chain(i-1,:);
    proposal = m_current + step_size .* randn(1, m_dim);

    % Check validity of the proposed model
    if any(isnan(doppler_fwd_chgdirection(proposal, Time_obs, 340))) || ...
       any(proposal < [0, 0, 100, 100, 0, 1, 10]) || ...
       any(proposal > [100, 300, 3000, 3000, 40, 179, 2000])
        chain(i,:) = m_current;
        continue
    end

    % Compute posterior probabilities
    log_post_current = log_posterior(m_current);
    log_post_proposal = log_posterior(proposal);

    % Metropolis-Hastings acceptance criterion
    alpha = exp(log_post_proposal - log_post_current);
    if rand < alpha
        chain(i,:) = proposal;
        accept_count = accept_count + 1;
    else
        chain(i,:) = m_current;
    end

    % Display acceptance rate every 100 iterations
    if mod(i, 100) == 0
        fprintf("Iter %d, acceptance = %.2f%%\n", i, 100 * accept_count / i);
    end
end

%% === Posterior Mean Prediction ===
m_mean = mean(chain(burn_in:end,:));
m_std = std(chain(burn_in:end,:));
f_pred = doppler_fwd_chgdirection(m_mean, Time_obs, 340);

% --- Plot observed vs predicted ---
figure;
plot(Time_obs, f_obs, 'ro', 'MarkerSize', 8); hold on
plot(Time_obs, f_pred, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
legend('Observed', 'Predicted (posterior mean)');
title('Observed vs Predicted Doppler Frequencies');

%% === Corner Plot for Posterior Samples ===
MC = chain(round(n_iter/2):end,:);  % Use second half of chain
MC = MC';  % Transpose to [n_param x N_samples]
param_names = {'f0', 'v0', 'L1', 'L2', 't0', 'theta', 'R'};
N = 7;
nbins = 60;
K = ones(3)/9;               % Smoothing kernel for 2D histograms
cmap = flipud(hot(150));     % Colormap

% Parameter plotting bounds
param_bounds = [ ...
    0, 100;    % f0
    1, 300;    % v0
    100, 5000; % L1
    100, 5000; % L2
    0, 50;     % t0
    1, 179;    % theta
    10, 1000   % R
];

for i = 1:length(m_mean)
    lower_bound = max(param_bounds(i,1), m_mean(i) - 3*m_std(i));
    upper_bound = min(param_bounds(i,2), m_mean(i) + 3*m_std(i));
    axis_limits{i} = [lower_bound, upper_bound];
end

% Precompute histogram bins
bins_all = cell(N,1);
bins_center = cell(N,1);
for i = 1:N
    bins = linspace(axis_limits{i}(1), axis_limits{i}(2), nbins);
    bins_all{i} = bins;
    bins_center{i} = (bins(1:end-1) + bins(2:end)) / 2;
end

% === Generate Corner Plot ===
figure('Color','w','Position',[100,100,1000,900]);
for i = 1:N
    for j = 1:i
        subplot(N,N,(i-1)*N + j)
        xi = MC(j,:);
        yi = MC(i,:);

        if i == j
            % Marginal distribution
            histogram(xi, bins_all{j}, 'Normalization','pdf', ...
                'FaceColor',[0.3,0.3,0.8], 'EdgeColor','none');
            xlim(axis_limits{j});
            title([param_names{i}, ' PDF']);
        else
            % Joint 2D distribution
            pdf = histcounts2(xi, yi, bins_all{j}, bins_all{i}, 'Normalization','pdf');
            pdf_smooth = conv2(pdf, K, 'same');
            pcolor(bins_center{j}, bins_center{i}, pdf_smooth'); shading flat
            colormap(gca, cmap);
            xlim(axis_limits{j});
            ylim(axis_limits{i});
        end

        % Axis labels
        if i == N, xlabel(param_names{j}); else, set(gca,'XTickLabel',[]); end
        if j == 1, ylabel(param_names{i}); else, set(gca,'YTickLabel',[]); end

        % Beautify
        set(gca, 'FontSize', 10, 'TickDir','out', 'Box','off');
    end
end
