%% 1. 数据准备 (Data Processing)
if ~exist('y', 'var') || ~exist('t', 'var')
    error('请先运行 main_include_T_2_b_new.m 生成数据 t 和 y');
end

% --- 数据插值 ---
data_new = zeros(100, 1001);
t_interp = 0:1:1000; 
for i = 1:100
    data_new(i,:) = interp1(t, y(:, i+1720), t_interp); 
end

% 聚合 Affinity 数据
data_kd = zeros(19, 1001);
for i = 1:100
    row_index = fix((i - 1)/10) - mod(i - 1, 10) + 10;
    for j = 1:1001
        data_kd(row_index, j) = data_kd(row_index, j) + data_new(i, j);
    end
end

% 网格定义
vals_kd = -26:1:-8;       
[TimeGrid, KdGrid] = meshgrid(t_interp, vals_kd); 

% --- 定义两种数据形式 ---
% 1. 线性数据 (用于 3D 几何形状)
FreqData_Linear = data_kd;

% 2. 对数数据 (用于 2D 热图颜色映射)
floor_val = 1e-6; % 噪音底限
FreqData_Log = data_kd;
FreqData_Log(FreqData_Log < floor_val) = floor_val;
FreqData_Log = log10(FreqData_Log);

% --- 提取病毒数据 ---
virus_idx = 2606; 
if size(y, 2) < virus_idx; virus_idx = 1822; end
virus_data = y(:, virus_idx);

%% 2. 绘图：多色系组合图 (Multi-Colormap Layout)
figure('Units', 'centimeters', 'Position', [4, 2, 18, 17], 'Color', 'w');
tlo = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ==========================================
% Panel a: Virus Load Dynamics (Line Plot)
% ==========================================
ax1 = nexttile([1 2]); % 横跨两列
p1 = plot(t, virus_data, 'LineWidth', 2.5);
p1.Color = [0.85, 0.33, 0.1]; % 砖红色

% 标题与标签
title(ax1, 'Viral Load Dynamics', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
% 添加 'a' 标签 (位于左上角外部)
text(ax1, -0.06, 1.1, 'a', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

ylabel(ax1, 'Virus Concentration', 'FontName', 'Arial', 'FontSize', 10);
xlabel(ax1, 'Time (hours)', 'FontName', 'Arial', 'FontSize', 10);
set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 9, 'FontName', 'Arial');
grid on; set(gca, 'GridAlpha', 0.15); xlim([0, max(t)]);

% ==========================================
% Panel b: 3D Clonal Landscape (Original Style)
% ==========================================
ax2 = nexttile;
% 使用线性数据绘制，保持物理高度
hSurf = surf(KdGrid, TimeGrid, FreqData_Linear);

% 样式设置
set(hSurf, 'EdgeColor', 'none', 'FaceColor', 'interp');
shading interp;
camlight('headlight'); lighting gouraud; material dull; 

% *** 关键点 1：为该坐标轴单独设置 colormap ***
colormap(ax2, parula); % 保持原来的 parula 色系

% 标题与标签
title(ax2, '3D Clonal Landscape', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
text(ax2, -0.15, 1.1, 'b', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

xlabel(ax2, 'Affinity (log K_d)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax2, 'Time', 'FontName', 'Arial', 'FontSize', 9);
zlabel(ax2, 'Frequency', 'FontName', 'Arial', 'FontSize', 8);

view(ax2, 35, 40); axis(ax2, 'tight');
set(ax2, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 9, 'FontName', 'Arial');

% 3D图的 Colorbar (线性)
cb2 = colorbar(ax2);
cb2.Label.String = 'Freq. (Linear)';
cb2.Label.FontSize = 8;
cb2.Layout.Tile = 'east'; % 让它紧贴子图

% ==========================================
% Panel c: Evolution Heatmap (Log Scale Color)
% ==========================================
ax3 = nexttile;
% 使用 Log 数据绘制颜色
hMap = pcolor(KdGrid, TimeGrid, FreqData_Log);
shading(ax3, 'interp'); set(hMap, 'EdgeColor', 'none');

% *** 关键点 2：为该坐标轴设置不同的 colormap ***
% 使用 turbo 或 jet 来增强 Log 数据的对比度
if exist('turbo', 'file')
    colormap(ax3, turbo); 
else
    colormap(ax3, jet); % 如果 Matlab 版本较老，回退到 jet
end

% 标题与标签
title(ax3, 'Affinity Maturation Heatmap', 'FontName', 'Arial', 'FontSize', 11, 'FontWeight', 'bold');
text(ax3, -0.15, 1.1, 'c', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 14, 'FontWeight', 'bold');

xlabel(ax3, 'Affinity (log K_d)', 'FontName', 'Arial', 'FontSize', 9);
ylabel(ax3, 'Time', 'FontName', 'Arial', 'FontSize', 9);

view(ax3, 2); axis(ax3, 'tight');
set(ax3, 'Layer', 'top', 'Box', 'on', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 9, 'FontName', 'Arial');

% Log Scale Colorbar 设置
cb3 = colorbar(ax3);
cb3.Label.String = 'Concentration (Log Scale)';
cb3.Label.FontSize = 8;

% 修正刻度显示为 10^x
c_limits = caxis(ax3); 
ticks = ceil(c_limits(1)):1:floor(c_limits(2));
cb3.Ticks = ticks;
cb3.TickLabels = arrayfun(@(x) sprintf('10^{%d}', fix(x)), ticks, 'UniformOutput', false);

%% 3. 导出建议
% print(gcf, 'Figure_Nature_Final.tif', '-dtiff', '-r600');