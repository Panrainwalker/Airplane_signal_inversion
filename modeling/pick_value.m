% -------------------------------------------------------------------------
% pick peak value 
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: Jun 17, 2025
% -------------------------------------------------------------------------
clc;clear;


%% Time-frequency
data = load(fullfile('FD2D/output', ['station', '3', '.txt']));
fs = 1 / (data(1,2) - data(1,1));  % 采样率
signal = data(2,:);         % 截取信号
dt = 1 / fs;
fc=5;


[s, f, t] = spectrogram(data(2,:), 256, 128, 256, fs, 'yaxis');
temp = log10(abs(s));
pcolor(t, f, temp);
shading interp;
clim([max(temp(:))/1.3, max(temp(:))*0.9])
colormap("turbo");
colorbar
xlabel('Time (s)')
ylabel('Frequency (Hz)')
set(gca,'FontSize',20)

ylim([fc*0.5 fc*2])
xlim([2 25])

disp('Click on points, press Enter when done.')
[clicked_x, clicked_y] = ginput;
plot(clicked_x, clicked_y, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
drawnow;

clicked_points = [clicked_x, clicked_y];
writematrix(clicked_points, 'clicked_points.txt', 'Delimiter', 'tab');
disp('Clicked points saved to clicked_points.txt')





