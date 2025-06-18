% -------------------------------------------------------------------------

% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: Jun 5, 2025
% -------------------------------------------------------------------------
%%%%%%

clc;clear;close all

%% parameters (Time, Stencil, grid, media, source)
%%%%% Time %%%%%
T_total = 30.0;     % unit (s)

%%%%% grid %%%%%
nz=500 + 6;
nx=1000 + 6;
length_z=2000;
length_x=4000; % length of x (m)

topo = 0 + zeros(1, nx); % topography

Flag_grid = 1; % 1: Cartesian Grid  /  2: Curvilinear Grid

Grid_disp = 0; 

%%%%% media %%%%%
Flag_media = 1;
Media_disp = 0;



%% Grid initial
Grid_init_t





%% Media initial
Media_init_t

CFL=1.2;
dt = ( CFL * dh / Vp_max) ;
% dt = 1e-5 ;
nt = round(T_total/dt) + 1;


%%

%%%%% source %%%%%
% -----------------------
% Moving source along x
% -----------------------
ndamp = 50;

fc = 5;                   % 固定频率 Hz
v0 = 120;                  % 源速度 m/s
src_iz = ndamp+30;  % 固定深度
src_ix_start = ndamp + 10;

% Source path in X
src_ix_start = ndamp+10;
src_ix_end = nx-ndamp+10;
src_ix_vec = src_ix_start : 1 : src_ix_end;  % All X positions the source moves through
Ns = length(src_ix_vec);  % Number of sources

t = (0:nt-1) * dt;
stf = sin(2 * pi * fc * t);   % 全时间长度的15Hz正弦波

%% Stencil 
Stencil_t;

RK4a = rk_coef{1}{1};
RK4b = rk_coef{1}{2};
RK6a = rk_coef{2}{1};
RK6b = rk_coef{2}{2};

half_fd_stentil = 3;
ni1 = half_fd_stentil+1;
ni2 = nx-half_fd_stentil;
nk1 = half_fd_stentil+1;
nk2 = nz-half_fd_stentil;

izvec = nk1 : nk2;     % interior z
ixvec = ni1 : ni2;     % interior x


MACF = [-0.30874,-0.6326,1.233,-0.3334,0.04168]; % forward
MACB = -MACF(end:-1:1); % backward

%% station 
stax = round(nx/2);
staz = nz-ndamp-300;
% staz = round(nz/2)
%% absorbing boundary
damp = ones(nz,nx);
ivec = ni1 : ni1+ndamp-1;
jvec = 1 : ndamp;
%-- z1/z2
for i = ni1 : ni2
    damp(ivec      ,i) = ((exp(-( (0.015*(ndamp-jvec)).^2 ) )).^10)';
    damp(nk2-ndamp+1:nk2,i) = ((exp(-( (0.015*(jvec      )).^2 ) )).^10)';
end
%-- x1/x2
for k = nk1 : nk2
    damp(k,ivec      ) = min( (exp(-( (0.015*(ndamp-jvec)).^2 ) )).^10, ...
                                 damp(k,ivec) );
    damp(k,ni2-ndamp+1:ni2) = min( (exp(-( (0.015*(jvec      )).^2 ) )).^10, ...
                                 damp(k,ni2-ndamp+1:ni2) );
end


%% solver

RK4a = RK4a .* dt;
RK4b = RK4b .* dt;
% RK6a = RK6a .* dt;
% RK6b = RK6b .* dt;


dx = xi_gd(2);
dz = abs(zt_gd(2));
MACF = MACF / dx;
MACB = MACB / dx;

CLT = CLT /dx;
sigma = (3*dt*Vs_min)/(4*dx); % for filter



flag_snap = 0;
gif_save= 1;
plot_sta = 1;
sta_name = '3';




% 运行 MAC 模型
sv_mac_rk_cart_ac_allstage;
Vzr_mac = Vzr;  % 保存结果

% sv_mac_rk_cart_allstage_noabs;
% Vzr_mac_noasb = Vzr;  % 保存结果


% sv_mac_rk_cart_os_allstage;
% Vzr_os = Vzr;  % 保存结果


% sv_ctl_rk_cart_allstage;
% sv_flt_rk_cart_allstage;






%%
% figure
% plot(Vzr)
if plot_sta ==1

figure
plot(Pr / max(abs(Pr)), 'r', 'DisplayName', 'MAC Model');
hold on
% plot(Vzr_mac_noasb / max(abs(Vzr_mac_noasb)), 'k--', 'DisplayName', 'OS Model');
% plot(10*(Vzr_os / max(abs(Vzr_os))-Vzr_mac / max(abs(Vzr_mac)))-1.5, 'b-', 'DisplayName', '10 * diff');

legend;
% ylim([-1, 1]);
xlabel('Time step');
ylabel('Normalized Vzr');
% title('Comparison of MAC and OS scheme');


fig = gcf;
fig.Units = 'inches';
fig.Position = [1, 1, 6, 4];



if ~exist('output', 'dir')
    mkdir('output');
end

cfl_str = strrep(sprintf('%.2f', CFL), '.', '_');


saveas(gcf, fullfile('output', ['station',sta_name, '.png']));
print(gcf, fullfile('output', ['station' ,sta_name,'.pdf']), '-dpdf', '-r300');

data = [t; Pr];
save(fullfile('output', ['station',sta_name, '.txt']), 'data', '-ascii');
end

%% Time-frequency
data = load(fullfile('output', ['station', sta_name, '.txt']));
fs = 1 / (data(1,2) - data(1,1));  % 采样率
signal = data(2,:);         % 截取信号
dt = 1 / fs;

% % === 补零 ===
% Npad = 2048;  % 可以根据需要设置为1024或更高
% signal_pad = [zeros(1,Npad), signal, zeros(1,Npad)];

% === 时间向量 ===
% tvec = linspace(0, length(signal_pad)*dt, length(signal_pad));

% === 绘图 ===
tiledlayout(2,1)

% --- 子图 1：补零后的原始信号 ---
nexttile
plot(data(1,:), data(2,:), 'k')
xlabel('Time (s)')
ylabel('Amplitude(count)')
title('Signal')
set(gca, 'FontSize', 20)
grid on
  xlim([2 25]);
% --- 子图 2：Spectrogram ---

% doppler effect
%four parameters 
L = (staz-src_iz)*(length_z/nz);

t0=14.75;
v0=v0;
fc=fc;
t=data(1,:);

%theory solution
t_ref = ( (t-t0) -sqrt( (t-t0).^2-(1-(v0/Vp(1,1)).^2).*((t-t0).^2-(L^2/Vp(1,1)^2))) ) ...
    ./ (1-(v0/Vp(1,1))^2);
f_sta = fc .* (1./(1+ ((v0./Vp(1,1)).*(v0.*t_ref./sqrt(L.^2+(v0.*t_ref).^2)) )));



nexttile
hold on
[s, f, t] = spectrogram(data(2,:), 256, 128, 256, fs, 'yaxis');
temp = log10(abs(s));
pcolor(t,f,temp);
h1 = plot(data(1,:), f_sta, 'w-', 'LineWidth', 4);  % 解析轨迹
h2 = yline(fc, 'w--', ['f_c = ', num2str(fc), ' Hz'], ...
    'LineWidth', 2, ...
    'FontSize', 18, ...
    'LabelHorizontalAlignment', 'left');
clim([max(temp(:))/1.3 max(temp(:))*0.9])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colormap("turbo");
colorbar
 ylim([fc*0.5 fc*2]);
  xlim([2 25]);
shading interp;
set(gca,'FontSize',20)
legend([h1 h2], { 'Analytical Solution'}, ...
    'Location', 'northeast', ...
    'FontSize', 16, ...
    'TextColor', 'w', ...
    'EdgeColor', 'k', ...           
    'Color', 'k', ...               
    'Box', 'on');  

annotation_text = {...
    ['$f_c$ = ', num2str(fc), ' Hz'], ...
    ['$v_0$ = ', num2str(v0), ' m/s'], ...
    ['$L$ = ', num2str(L), ' m'], ...
    ['$t_0$ = ', num2str(t0), ' s']};

annotation('textbox', [0.8, 0.3, 0.1, 0.1], ...
    'String', annotation_text, ...
    'EdgeColor', 'w', ...               % 白色边框
    'LineWidth', 1.5, ...               % 可选：加粗边框
    'FontSize', 16, ...
    'Color', 'k', ...
    'BackgroundColor', 'w', ...
    'Interpreter', 'latex');

























function cmap = seismic_colormap(m)
    if nargin < 1
        m = 256; % 默认256色
    end
    % 线性蓝->白->红分三段生成
    n = floor(m/2);
    blue = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)];    % 蓝到白
    red  = [ones(m-n,1), linspace(1,0,m-n)', linspace(1,0,m-n)']; % 白到红
    cmap = [blue; red];
end