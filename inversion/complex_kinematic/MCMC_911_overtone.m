clc; clear; close all;
addpath(fullfile(pwd, '../modeling'));

%% 1. Load four observed 911 overtones
data1 = load('../../data/911_overtone1.txt');  
data2 = load('../../data/911_overtone2.txt');  
data3 = load('../../data/911_overtone3.txt');  
data4 = load('../../data/911_overtone4.txt');  

Time1 = data1(:,1); f1 = data1(:,2);
Time2 = data2(:,1); f2 = data2(:,2);
Time3 = data3(:,1); f3 = data3(:,2);
Time4 = data4(:,1); f4 = data4(:,2);

%% 2. Prior settings
prior_mu  = [80.46 214.92 357.78  574.25 193.80 411.47 4.84 5.36];
prior_std = [1, 1, 1, 1, 1, 1, 1, 1];  % Can be tuned empirically
log_prior = @(m) -0.5 * sum(((m - prior_mu)./prior_std).^2);

%% 3. MCMC setup
n_iter = 200000;
burn_in = n_iter * 0.4;
m_dim = 8;  % Parameters: f0_1, f0_2, f0_3, f0_4, v0, L, t0, a
chain = zeros(n_iter, m_dim);

% Initialization
chain(1,:) = prior_mu;  
step_size = [5,5,5,5,2,2,1,1]/200;

%% 4. Log-posterior definition
log_likelihood = @(m) ...
    -0.5 * sum((f1 - doppler_fwd_acc([m(1), m(5:8)], Time1, 340)).^2,'omitnan') ...
    -0.5 * sum((f2 - doppler_fwd_acc([m(2), m(5:8)], Time2, 340)).^2,'omitnan') ...
    -0.5 * sum((f3 - doppler_fwd_acc([m(3), m(5:8)], Time3, 340)).^2,'omitnan') ...
    -0.5 * sum((f4 - doppler_fwd_acc([m(4), m(5:8)], Time4, 340)).^2,'omitnan');

log_posterior = @(m) log_prior(m) + log_likelihood(m);

%% 5. MCMC iteration loop
accept_count = 0;
for i = 2:n_iter
    m_curr = chain(i-1,:);
    proposal = m_curr + step_size .* randn(1,m_dim);

    % Validity check
    if any(proposal(1:4) < 0 | proposal(1:4) > 1000) || ...
       any(proposal(5:8) < [1,1,1,-3]) || any(proposal(5:8) > [400,1000,10,10])
        chain(i,:) = m_curr;
        continue
    end

    log_post_curr = log_posterior(m_curr);
    log_post_prop = log_posterior(proposal);

    alpha = exp(log_post_prop - log_post_curr);
    if rand < alpha
        chain(i,:) = proposal;
        accept_count = accept_count + 1;
    else
        chain(i,:) = m_curr;
    end

    if mod(i,200) == 0
        fprintf('Iter %d, acceptance = %.2f%%\n', i, 100*accept_count/i);
    end
end

%% 6. Posterior analysis
post = chain(round(burn_in):end,:);
f0_post = post(:,1:4);
shared_post = post(:,5:8);

fprintf('Posterior mean f0_1..4 = %.2f %.2f %.2f %.2f\n', mean(f0_post));
fprintf('Shared parameters mean: v0=%.2f, L=%.2f, t0=%.2f, a=%.2f\n', mean(shared_post));

%% 7. Posterior mean prediction
m_mean = mean(chain(burn_in:end,:));
m_std = std(chain(burn_in:end,:), 0, 1);  % Standard deviation of each parameter
f_pred1 = doppler_fwd_acc([m_mean(1), m_mean(5:8)], Time1, 340);
f_pred2 = doppler_fwd_acc([m_mean(2), m_mean(5:8)], Time2, 340);
f_pred3 = doppler_fwd_acc([m_mean(3), m_mean(5:8)], Time3, 340);
f_pred4 = doppler_fwd_acc([m_mean(4), m_mean(5:8)], Time4, 340);

%% 8. Plot observed and predicted frequency-time curves
figure; hold on;

Time = m_mean(7)-4.5:0.01:m_mean(7)+4;
f_pred1 = doppler_fwd_acc([m_mean(1), m_mean(5:8)], Time, 340);
f_pred2 = doppler_fwd_acc([m_mean(2), m_mean(5:8)], Time, 340);
f_pred3 = doppler_fwd_acc([m_mean(3), m_mean(5:8)], Time, 340);
f_pred4 = doppler_fwd_acc([m_mean(4), m_mean(5:8)], Time, 340);

% Observations
plot(Time1, f1, 'ro','MarkerSize',6,'DisplayName','Observed 1');
plot(Time2, f2, 'go','MarkerSize',6,'DisplayName','Observed 2');
plot(Time3, f3, 'bo','MarkerSize',6,'DisplayName','Observed 3');
plot(Time4, f4, 'mo','MarkerSize',6,'DisplayName','Observed 4');

% Posterior mean predictions
plot(Time, f_pred1, 'r-','LineWidth',2,'DisplayName','Predicted 1');
plot(Time, f_pred2, 'g-','LineWidth',2,'DisplayName','Predicted 2');
plot(Time, f_pred3, 'b-','LineWidth',2,'DisplayName','Predicted 3');
plot(Time, f_pred4, 'm-','LineWidth',2,'DisplayName','Predicted 4');

xlabel('Time (s)');
ylabel('Frequency (Hz)');
legend('Location','best');
title('Observed and Predicted Frequencies for Four Data Sets');
grid on;
set(gca,'YScale','log');

%% 9. Save predicted and observed data
% save('911_pred1.txt', 'f_pred1', '-ascii');
% save('911_pred2.txt', 'f_pred2', '-ascii');
% save('911_pred3.txt', 'f_pred3', '-ascii');
% save('911_pred4.txt', 'f_pred4', '-ascii');
% 
% save('911_obs1.txt', 'f1', '-ascii');
% save('911_obs2.txt', 'f2', '-ascii');
% save('911_obs3.txt', 'f3', '-ascii');
% save('911_obs4.txt', 'f4', '-ascii');

%% 10. Corner plot for 8D MCMC posterior samples
MC = chain(round(n_iter/2):end,:);  % [N_samples x n_param]
MC = MC';  % Transpose to [n_param x N_samples]

param_names = {'f0_1','f0_2','f0_3','f0_4','v0','L','t0','a'};
N = length(param_names);
nbins = 60;
K = ones(3)/9;        % Smoothing kernel
cmap = flipud(hot(150));  % Color map

% Parameter bounds (adjustable)
param_bounds = [ ...
    0, 1000;    % f0_1
    0, 1000;    % f0_2
    0, 1000;    % f0_3
    0, 1000;    % f0_4
    1, 400;     % v0
    1, 1000;    % L
    0, 10;      % t0
    -3, 10      % a
];

% Compute axis limits
for i = 1:N
    lower_bound = max(param_bounds(i,1), m_mean(i) - 3*m_std(i));
    upper_bound = min(param_bounds(i,2), m_mean(i) + 3*m_std(i));
    axis_limits{i} = [lower_bound, upper_bound];
end

% Precompute histogram bins
bins_all = cell(N,1); bins_center = cell(N,1);
for i = 1:N
    bins = linspace(axis_limits{i}(1), axis_limits{i}(2), nbins);
    bins_all{i} = bins;
    bins_center{i} = (bins(1:end-1) + bins(2:end))/2;
end

% Generate corner plot
figure('Color','w','Position',[100,100,1200,1200]);
for i = 1:N
    for j = 1:i
        subplot(N,N,(i-1)*N + j)
        xi = MC(j,:);
        yi = MC(i,:);

        if i == j
            % Marginal PDFs
            histogram(xi, bins_all{j}, 'Normalization','pdf', ...
                'FaceColor',[0.3,0.3,0.8], 'EdgeColor','none');
            xlim(axis_limits{j});
            title([param_names{i}, ' PDF']);
        else
            % Joint PDFs
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

        % Aesthetic adjustments
        set(gca, 'FontSize', 10, 'TickDir','out', 'Box','off');
    end
end
