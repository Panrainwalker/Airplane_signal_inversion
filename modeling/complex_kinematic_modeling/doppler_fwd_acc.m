% -------------------------------------------------------------------------
% doppler effect analytical solution (forward)
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: Jun 17, 2025
% -------------------------------------------------------------------------

function [f_sta] = doppler_fwd_acc(m_pred,t,Vp)

fc = m_pred(1);            % Source frequency (Hz)
v0 = m_pred(2);          % Horizontal rupture velocity (m/s)
L = m_pred(3);           % Perpendicular distance from source to station (m)
t0 = m_pred(4);        % Source activation time (s)
a = m_pred(5); 


% % Piecewise source position and velocity
% x_func = @(t_emit)  ( v0 * (t_emit-t0) + 0.5 * a * (t_emit-t0).^2) ;
% % x_func = @(t_emit)  (x0 + v0 * t_emit + 0.5 * a * t_emit.^2);
% 
% v_func = @(t_emit)  (v0 + a * (t_emit-t0)) ;
% 
% x_func = @(t_emit) (t_emit >= t0) .* (v0 * (t_emit - t0) + 0.5 * a * (t_emit - t0).^2);
% v_func = @(t_emit) (t_emit >= t0) .* (v0 + a * (t_emit - t0));
% % v_func = @(t_emit)  (v0 + a * t_emit) ;
% % R function: distance from source to station
% R_func = @(t_emit) sqrt((x_func((t_emit))).^2 + L^2);


x_func = @(t_emit)  ( v0 * (t_emit-t0) + 0.5 * a * (t_emit-t0).^2) ;
v_func = @(t_emit)  (v0 + a * (t_emit-t0)) ;
R_func = @(t_emit) sqrt(x_func(t_emit).^2 + L^2);

% Invert travel time equation to find emission time
t_emit = NaN(size(t));
for i = 1:length(t)
    t_obs = t(i);
    travel_eq = @(t_r) t_obs - t_r - R_func(t_r)/Vp;
    try
        t_emit(i) = fzero(travel_eq, [t0-100, t_obs]);
    catch
        t_emit(i) = NaN;
    end
end

% Evaluate functions
x_val = x_func(t_emit);
v_val = v_func(t_emit);
v_radial = v_val .* x_val ./ sqrt(x_val.^2 + L^2);

f_sta = fc ./ (1 + v_radial / Vp);  % Doppler shift formula

end




