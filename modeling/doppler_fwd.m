% -------------------------------------------------------------------------
% doppler effect analytical solution (forward)
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: Jun 17, 2025
% -------------------------------------------------------------------------

function [f_sta] = doppler_fwd(m_pred,t,Vp)

fc = m_pred(1);            % Source frequency (Hz)
v0 = m_pred(2);          % Horizontal rupture velocity (m/s)
L = m_pred(3);           % Perpendicular distance from source to station (m)
t0 = m_pred(4);        % Source activation time (s)


t_ref = ((t - t0) - sqrt((t - t0).^2 - (1 - (v0 / Vp)^2) .* ((t - t0).^2 - (L^2 / Vp^2)))) ...
        ./ (1 - (v0 / Vp)^2);


f_sta = fc .* (1 ./ (1 + ((v0 / Vp) .* (v0 .* t_ref ./ sqrt(L^2 + (v0 .* t_ref).^2)))));


end



