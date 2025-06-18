% -------------------------------------------------------------------------
% doppler effect analytical solution
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: Jun 17, 2025
% -------------------------------------------------------------------------



clc; clear;

% Doppler Effect Simulation for a Moving  Source
% Define source parameters
L = 300;           % Perpendicular distance from source to station (m)
t0 = 15.75;        % Source activation time (s)
v0 = 120;          % Horizontal rupture velocity (m/s)
fc = 5;            % Source frequency (Hz)
Vp = 340;          % Propagation velocity (m/s)
t = 0:0.001:30;   % Time vector (s)

% Analytical solution for the observed frequency at the station
% Reference time accounting for source motion and geometry
t_ref = ((t - t0) - sqrt((t - t0).^2 - (1 - (v0 / Vp)^2) .* ((t - t0).^2 - (L^2 / Vp^2)))) ...
        ./ (1 - (v0 / Vp)^2);

% Observed frequency at the station as a function of time
f_sta = fc .* (1 ./ (1 + ((v0 / Vp) .* (v0 .* t_ref ./ sqrt(L^2 + (v0 .* t_ref).^2)))));

% Plot the analytical solution
figure;
plot(t, f_sta, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Observed Frequency (Hz)');
grid on;
set(gca, 'FontSize', 14);