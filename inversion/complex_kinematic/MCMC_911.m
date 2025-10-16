clc; clear; close all;

% --- Load observed data ---
data = load('../../data/911.txt');
% data = load('../modeling/clicked_points.txt');
Time_obs = data(:,1);
f_obs = data(:,2);

w_obs = ones(size(f_obs));  % Default weight for each observation = 1

% --- Prior distribution parameters ---
prior_mu = [223.5473, 177.2714, 485.3236, 5.0034, 1.8430];
prior_std = [200, 100, 500, 2, 3];

% log_prior = @(m) sum(log(normpdf(m, prior_mu, prior_std)));
log_prior = @(m) -0.5 * sum(((m - prior_mu)./prior_std).^2);

% --- Likelihood function ---
log_likelihood = @(m) ...
    -0.5 * sum(w_obs .* (f_obs - doppler_fwd_acc(m, Time_obs, 340)).^2, 'omitnan');

% --- Posterior function ---
log_posterior = @(m) log_prior(m) + log_likelihood(m);

% --- MCMC settings ---
n_iter = 20000;
burn_in = n_iter * 0.4;
m_dim = 5;
chain = zeros(n_iter, m_dim);
chain(1,:) = prior_mu;  % Initialize chain

step_size = [10, 10, 50, 1, 2] / 200;  % Initial step size

accept_count = 0;
adapt_interval = 500;   % Interval for potential adaptive tuning

% --- Metropolis-Hastings sampling loop ---
for i = 2:n_iter
    m_current = chain(i-1,:);
    proposal = m_current + step_size .* randn(1, m_dim);

    % Validity check: avoid invalid or NaN model outputs
    if any(isnan(doppler_fwd_acc(proposal, Time_obs, 340))) || ...
       any(proposal < [10, 1, 10, 1, -10]) || any(proposal > [400, 1000, 2000, 60, 10])
        chain(i,:) = m_current;
        continue
    end

    % Compute posterior probabilities
    log_post_current = log_posterior(m_current);
    log_post_proposal = log_posterior(proposal);

    % Metropolis acceptance criterion
    alpha = exp(log_post_proposal - log_post_current);
    if rand < alpha
        chain(i,:) = proposal;
        accept_count = accept_count + 1;
    else
        chain(i,:) = m_current;
    end

    % Print progress every 100 iterations
    if mod(i, 100) == 0
        acceptance_rate = accept_count / i;
        fprintf("Iter %d, acceptance = %.2f%%\n", i, 100 * acceptance_rate);
    end
end

%% --- Posterior mean prediction ---
m_mean = mean(chain(burn_in:end,:));
m_std = std(chain(burn_in:end,:));
f_pred = doppler_fwd_acc(m_mean, Time_obs, 340);

% --- Plot observed vs. predicted frequencies ---
figure;
plot(Time_obs, f_obs, 'ro', 'MarkerSize', 8); hold on
plot(Time_obs, f_pred, 'k-', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Frequency (Hz)');
legend('Observed', 'Predicted (posterior mean)');

%% --- Corner plot for MCMC posterior samples ---
MC = chain(round(n_iter/2):end,:);  % [N_samples x n_param]
MC = MC';  % Transpose to [n_param x N_samples]
param_names = {'f0', 'v0', 'L', 't0', 'a'};
N = 5;
nbins = 60;
K = ones(3)/9;        % Smoothing kernel
cmap = flipud(hot(150));  % Color map

param_bounds = [ ...
    0, 5000;    % f0
    0, 2000;    % v0
    1, 5000;    % L
    0, 10000;   % t0
    -10, 10     % a
];

% Determine axis limits based on posterior means and stds
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

% --- Draw corner plot ---
figure('Color','w','Position',[100,100,1000,900]);
for i = 1:N
    for j = 1:i
        subplot(N, N, (i-1)*N + j)
        xi = MC(j,:);
        yi = MC(i,:);

        if i == j
            % Marginal posterior PDF
            histogram(xi, bins_all{j}, 'Normalization','pdf', ...
                'FaceColor',[0.3,0.3,0.8], 'EdgeColor','none');
            xlim(axis_limits{j});
            title([param_names{i}, ' PDF']);
        else
            % Joint posterior PDF
            pdf = histcounts2(xi, yi, bins_all{j}, bins_all{i}, 'Normalization','pdf');
            pdf_smooth = conv2(pdf, K, 'same');
            pcolor(bins_center{j}, bins_center{i}, pdf_smooth'); shading flat
            colormap(gca, cmap);
            xlim(axis_limits{j});
            ylim(axis_limits{i});
        end

        % Axis labeling
        if i == N, xlabel(param_names{j}); else, set(gca,'XTickLabel',[]); end
        if j == 1, ylabel(param_names{i}); else, set(gca,'YTickLabel',[]); end

        % Beautify
        set(gca, 'FontSize', 10, 'TickDir','out', 'Box','off');
    end
end
