% -------------------------------------------------------------------------
% solver of acoustic using cartesian grid
% all stage forward calculation for right hand side (rhs)
%
%
% Author: Yuxing Pan (panyx2023@mail.sustech.edu.cn)
% Affiliation: Southern University of Science and Technology (SUSTech)
% Date: May 31, 2025
% -------------------------------------------------------------------------
%%%%%%

%% predefine
it = 0;

hWz = zeros(nz,nx,3); % rhs : [Vx_z Vz_z P_z]
hWx = zeros(nz,nx,3); % rhs : [Vx_x Vz_x P_x]
hW = zeros(nz,nx,3);  % rhs : [Vx_t Vz_t P_t]
W = zeros(nz,nx,3);   % rhs : [Vx Vz P]
tW = zeros(nz,nx,3);   % middle variance
mW = zeros(nz,nx,3);   % middle variance


while 1

    it = it + 1
    if it > nt
        break
    end
    mW = W;

    % Stage 1
    W = mW;
    Add_source;
    apply_abs_ac
    % prep_free_cart;

    % for inner point
    cal_macF
    % for free surface
    % apply_free_srf_cart
    ass_rhs_cart_ac;


    k1 = hW;


    % Stage 2
    W = mW + RK4a(2) * k1;
    Add_source;
    apply_abs_ac
    %prep_free_cart;
    % for inner point
    cal_macB
    % for free surface
    %apply_free_srf_cart
    ass_rhs_cart_ac;

    k2 = hW;


    % Stage 3
    W = mW + RK4a(3) * k2;
    Add_source;
    apply_abs_ac
    %prep_free_cart;
    % for inner point
    cal_macF
    % for free surface
    %apply_free_srf_cart
    ass_rhs_cart_ac;


    k3 = hW;



    % Stage 4
    W = mW + RK4a(4) * k3;
    Add_source;
    apply_abs_ac
    %prep_free_cart;
    % for inner point
    cal_macB
    % for free surface
    %apply_free_srf_cart
    ass_rhs_cart_ac;


    k4 = hW;


    % Final RK4 combination
    W = mW + RK4b(1) * k1 + RK4b(2) * k2 + RK4b(3) * k3 + RK4b(4) * k4;
    apply_abs_ac

    Vxr(it) = W(staz  ,stax,1);
    Vzr(it) = W(staz  ,stax,2);
    Pr(it) = W(staz  ,stax,3);
    it = it + 1
    if it > nt
        break
    end
    mW = W;


    % Stage 1
    W = mW;
    Add_source;
    apply_abs_ac
    %prep_free_cart;
    % for inner point
    cal_macB
    % for free surface
    %apply_free_srf_cart
    ass_rhs_cart_ac;

    % for free surface

    k1 = hW;


    % Stage 2
    W = mW + RK4a(2) * k1;
    Add_source;
    apply_abs_ac
    %prep_free_cart;
    % for inner point
    cal_macF
    % for free surface
    %apply_free_srf_cart
    ass_rhs_cart_ac;


    k2 = hW;


    % Stage 3
    W = mW + RK4a(3) * k2;
    Add_source;
    apply_abs_ac
    %prep_free_cart;
    % for inner point
    cal_macB
    % for free surface
    %apply_free_srf_cart
    ass_rhs_cart_ac;


    k3 = hW;


    % Stage 4
    W = mW + RK4a(4) * k3;
    Add_source;
    apply_abs_ac
    %prep_free_cart;
    % for inner point
    cal_macF
    % for free surface
    %apply_free_srf_cart
    ass_rhs_cart_ac;


    k4 = hW;


    % Final RK4 combination
    W = mW + RK4b(1) * k1 + RK4b(2) * k2 + RK4b(3) * k3 + RK4b(4) * k4;
    apply_abs_ac


    Vxr(it) = W(staz  ,stax,1);
    Vzr(it) = W(staz  ,stax,2);
    Pr(it) = W(staz  ,stax,3);


if flag_snap == 1 && mod(it, 100) == 0
    figure(99); clf;
    pcolor(x_gd, z_gd, W(:,:,3));
    shading interp;
colormap(seismic_colormap(256));
colorbar;

    set(gca, 'YDir', 'reverse');
    axis equal tight;
    caxis([-5e7 5e7]);

    title(['P: ', num2str(it * dt, '%.3f'), ' s']);
    xlabel('X (m)');
    ylabel('Z (m)');
    hold on;

    % 台站位置
    h2 = plot(x_gd(1, stax), z_gd(staz, 1), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');

    % 当前源位置
    ix = round(src_ix_start + v0 * it * dt / dx);
    ix = min(max(ix, 1), size(x_gd, 2));
    h3 = plot(x_gd(1, ix), z_gd(src_iz, 1), 'k*', 'MarkerSize', 10,'MarkerFaceColor', 'y');

  
    h1 = plot(nan, nan, 'k');
legend([h2, h3], {'Station', 'Source'}, ...
       'Location', 'northwest', ...
       'FontSize', 12, ...
       'Box', 'on');  % 可加 Box 显示边框
hold off;
    drawnow;
end

if gif_save == 1 && mod(it, 100) == 0
    % 创建输出文件夹
    png_dir = fullfile('output', 'frames');
    if ~exist(png_dir, 'dir')
        mkdir(png_dir);
    end

    % 创建不可见图窗
    fig = figure('Visible', 'off', 'Position', [100, 100, 800, 600]);

    % 绘图
    ax = axes('Parent', fig);
    pcolor(ax, x_gd, z_gd, W(:,:,3));
    shading(ax, 'interp');
    colormap(ax, seismic_colormap(256));
    caxis(ax, [-5e7 5e7]);
    set(ax, 'YDir', 'reverse');
    axis(ax, 'equal', 'tight');
    title(ax, ['P: ', num2str(it * dt, '%.3f'), ' s']);
    xlabel(ax, 'X (m)');
    ylabel(ax, 'Z (m)');
    hold(ax, 'on');

    % 台站和震源
    plot(ax, x_gd(1, stax), z_gd(staz, 1), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
    ix = round(src_ix_start + v0 * it * dt / dx);
    ix = min(max(ix, 1), size(x_gd, 2));
    plot(ax, x_gd(1, ix), z_gd(src_iz, 1), 'k*', 'MarkerSize', 10, 'MarkerFaceColor', 'y');

        % 台站位置
    h2 = plot(x_gd(1, stax), z_gd(staz, 1), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
    h3 = plot(x_gd(1, ix), z_gd(src_iz, 1), 'k*', 'MarkerSize', 10,'MarkerFaceColor', 'y');

    legend([h2, h3], {'Station', 'Source'}, ...
       'Location', 'northwest', ...
       'FontSize', 12, ...
       'Box', 'on');  % 可加 Box 显示边框
hold off;
    % 保存 PNG
    png_name = sprintf('frame_%05d.png', it);
    png_path = fullfile(png_dir, png_name);
    exportgraphics(fig, png_path, 'Resolution', 150);

    % 关闭图窗
    close(fig);
end
end

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