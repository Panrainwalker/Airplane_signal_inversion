% -------------------------------------------------------------------------
% apply IRLS to parameters estimation
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: Jun 18, 2025
% -------------------------------------------------------------------------

clc;clear;
addpath(fullfile(pwd, '../modeling'));

%% initial guess
m = zeros(1,4);
m(1) = 6 ; % freq
m(2) = 100 ; % source vel
m(3) = 400 ; % dist
m(4) = 20 ; %  t0

t = -50:0.001:50;

% forward for initial guess
f_pred = doppler_fwd(m,t,340);

%% observation (use FD result)
data = load('../modeling/clicked_points.txt');
Time_obs =  data(:,1);
f_obs = data(:,2);      

%% apply IRLS
m_mid_f = zeros(size(m));   % to get derivatives
m_mid_b = zeros(size(m));   % to get derivatives
m_updata = zeros(size(m));% updata models


% compare inital guess with obs

f_pred_init = doppler_fwd(m,Time_obs,340);


nit=10000; %number of iteration

% start IRLS
for it=1:nit

f_pred = doppler_fwd(m,Time_obs,340);
diff = f_obs - f_pred;    % difference between obs and pred

Gm = zeros(length(Time_obs),1); % G for four models (nt,1)
alpha = 0.01 ; % increment
loss_record(it) = sum(diff.^2);
m_mid_f = m;
m_mid_b = m;
for i = 1:4
    m_mid_f(i)=m(i)+alpha*m(i);
    m_mid_b(i)=m(i)-alpha*m(i);
    Gm = (doppler_fwd(m_mid_f,Time_obs,340) - doppler_fwd(m_mid_b,Time_obs,340))...
                    ./ (2.*alpha.*m(i));
    m_updata(i) = m(i)' + inv(Gm'*Gm)*Gm'*(diff);
end


m = m_updata';

end

f_pred = doppler_fwd(m, Time_obs, 340);

%%

figure
hold on
plot(Time_obs, f_obs, 'ro', 'MarkerSize', 8, 'LineWidth', 2)
plot(Time_obs, f_pred, 'k-', 'LineWidth', 2)
plot(Time_obs, f_pred_init, 'k--', 'LineWidth', 2)

legend({'Observed', 'Fitted', 'Initial Guess'}, ...
    'FontSize', 16, ...
    'Location', 'northeast');

xlabel('Time (s)', 'FontSize', 18)
ylabel('Frequency (Hz)', 'FontSize', 18)
set(gca, 'FontSize', 20)

% 注释文本内容（LaTeX 样式更美观）
annotation_text = {...
    ['Predicted parameters'], ...
    [''], ...
    ['$f_c$ = ', num2str(m_updata(1)), ' Hz'], ...
    ['$v_0$ = ', num2str(m_updata(2)), ' m/s'], ...
    ['$L$ = ', num2str(m_updata(3)), ' m'], ...
    ['$t_0$ = ', num2str(m_updata(4)), ' s']};

% 添加注释框
annotation('textbox', [0.15, 0.3, 0.15, 0.2], ...
    'String', annotation_text, ...
    'EdgeColor', 'k', ...        
    'LineWidth', 1.2, ...
    'FontSize', 20, ...
    'Color', 'k', ...
    'BackgroundColor', 'w', ...
    'Interpreter', 'latex');      


