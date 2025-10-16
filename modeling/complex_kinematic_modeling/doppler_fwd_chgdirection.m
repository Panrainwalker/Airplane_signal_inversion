function f_sta_at_obs_t = doppler_fwd_chgdirection_new(m_pred, t_obs, c)
% -------------------------------------------------------------------------
% doppler effect analytical solution (forward) - 无 x_start 依赖版本（首段自动延长）
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn), 改进 by ChatGPT
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: Jun 17, 2025 (Updated: Aug 10, 2025)
% -------------------------------------------------------------------------
% 输入：
%   m_pred: 参数向量 [f0, v_mag, L1, L2, t0, theta_end_deg, R]
%   t_obs:  观测时间点，列向量
%   c:      波速 (m/s)
% 输出：
%   f_sta_at_obs_t: 观测时间点对应的多普勒频率，列向量

f0 = m_pred(1);             % Source frequency (Hz)
v_mag = m_pred(2);          % Horizontal rupture velocity (m/s)
L1 = m_pred(3);             % 垂直距离段1
L2 = m_pred(4);             % 垂直距离段2
t0 = m_pred(5);             % Source activation time (s)
theta_end_deg = m_pred(6);  % 转弯角度 (deg)
R = m_pred(7);              % Arc radius

%% 基础参数
N_seg = 200;
theta_end = deg2rad(theta_end_deg);
P0 = [0; 0];           % 接收器位置（原点）
ref_point = [0; L1];   % 轨迹第一段的参考点

% 圆心坐标
C_x = (L1 + (L1 - L2) / (cos(theta_end) - 1)) * tan(theta_end/2) - R * tan(theta_end/2);
C = P0 + [C_x; L1 - R];

% 圆弧轨迹
theta_arc = linspace(deg2rad(90), deg2rad(90 - theta_end_deg), N_seg);
arc = R * [cos(theta_arc); sin(theta_arc)] + C;
arc_start = arc(:, 1);

%% 第一段直线（自动延长到提前 50 s）
T_extend = 50; % s
L_extend = v_mag * T_extend;
dir_line1 = (arc_start - ref_point);
dir_line1 = dir_line1 / norm(dir_line1);
start_point = ref_point - dir_line1 * L_extend;
line1 = [linspace(start_point(1), arc_start(1), N_seg);
         linspace(start_point(2), arc_start(2), N_seg)];

%% 第三段直线（沿末端切线）
theta_end_rad = theta_arc(end);
tangent_dir = [-sin(theta_end_rad); cos(theta_end_rad)];
tangent_start = arc(:, end);
L3 = 2000; % 延伸长度
line3 = [linspace(tangent_start(1), tangent_start(1) - tangent_dir(1)*L3, N_seg);
         linspace(tangent_start(2), tangent_start(2) - tangent_dir(2)*L3, N_seg)];

%% 合并轨迹（不再用 x_start 平移）
trajectory = [line1, arc, line3];

%% 计算速度向量
vel = diff(trajectory, 1, 2);
vel = [vel, vel(:, end)];
t_hat = vel ./ vecnorm(vel);
v = v_mag * t_hat;

%% 计算接收时间
dist_to_receiver = vecnorm(trajectory - P0, 2, 1);
ds = vecnorm(diff(trajectory, 1, 2));
dt = ds / v_mag;
t_emission = [0, cumsum(dt)];
t_receive_raw = t_emission + dist_to_receiver / c;

% 排序 + 去重
[t_receive_sorted, sort_idx] = sort(t_receive_raw);
trajectory_sorted = trajectory(:, sort_idx);
v_sorted = v(:, sort_idx);
[t_receive_unique, unique_idx] = unique(t_receive_sorted, 'stable');
trajectory_unique = trajectory_sorted(:, unique_idx);
v_unique = v_sorted(:, unique_idx);

%% 自动找 (0, L1) 作为时间参考点
[~, idx_ref] = min(vecnorm(trajectory_unique - ref_point, 2, 1));
t_receive_ref = t_receive_unique(idx_ref);
t_receive_shifted = t_receive_unique - t_receive_ref;

%% 插值（以 t0 为参考）
traj_at_obs_t = trajectory_at_time(t_obs - t0, trajectory_unique, t_receive_shifted);
v_at_obs_t    = interpolate_velocity_at_time(t_obs - t0, t_receive_shifted, v_unique);

% 填补 NaN
for dim = 1:2
    traj_at_obs_t(dim, :) = fillmissing(traj_at_obs_t(dim, :), 'linear', 'EndValues', 'extrap');
    v_at_obs_t(dim, :) = fillmissing(v_at_obs_t(dim, :), 'linear', 'EndValues', 'extrap');
end

%% 计算法向速度
v_normal_obs = zeros(length(t_obs), 1);
for i = 1:length(t_obs)
    traj_pt = traj_at_obs_t(:, i);
    v_pt = v_at_obs_t(:, i);
    n = traj_pt - P0;
    n_hat = n / norm(n);
    v_normal_obs(i) = dot(v_pt, n_hat);
end

%% Doppler 频率
f_sta_at_obs_t = f0 ./ (1 - (-v_normal_obs / c));

end

%% 子函数
function traj_at_t = trajectory_at_time(t_obs_shifted, trajectory, t_receive_sorted)
    t_obs_shifted = t_obs_shifted(:);
    traj_x = interp1(t_receive_sorted, trajectory(1,:), t_obs_shifted, 'linear', 'extrap');
    traj_y = interp1(t_receive_sorted, trajectory(2,:), t_obs_shifted, 'linear', 'extrap');
    traj_at_t = [traj_x.'; traj_y.'];
end

function v_at_t = interpolate_velocity_at_time(t_obs_shifted, t_receive_sorted, v)
    t_obs_shifted = t_obs_shifted(:);
    vx = interp1(t_receive_sorted, v(1,:), t_obs_shifted, 'linear', 'extrap');
    vy = interp1(t_receive_sorted, v(2,:), t_obs_shifted, 'linear', 'extrap');
    v_at_t = [vx.'; vy.'];
end
