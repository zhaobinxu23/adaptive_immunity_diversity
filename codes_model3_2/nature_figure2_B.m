clc; clear; close all;

%% ==========================
% 1. 模型参数与模拟计算
% ==========================

% --- 模拟 A: 对照组 (Normal People) ---
x0_norm(1) = 1e5;   % Antibody
x0_norm(2) = 1e10;  % Vaccine Antigen
x0_norm(3) = 0;     % Complex
x0_norm(4) = 1000;  % Tc cells

para_norm(1) = 1e-7; 
para_norm(2) = 1e-14;
para_norm(3) = 5;      % [Normal] Feedback constant
para_norm(4) = 2; 
para_norm(5) = -0.05; 
para_norm(6) = 1e-5; 
para_norm(7) = 1e3;    
para_norm(8) = 0.01;
para_norm(9) = 10;
para_norm(10)= 0.01;

% 求解器设置
options = odeset('RelTol',1e-5, 'AbsTol',1e-6);
[t_norm, y_norm] = ode15s(@pathway_model_4_antibody_immune_cellular_humoral, ...
                          [0 500], x0_norm, options, para_norm);

% --- 模拟 B: 实验组 (Rituximab / Blocked Humoral) ---
x0_exp = x0_norm;

para_exp = para_norm;
para_exp(3) = 3;  % [Treatment] 降低抗体反馈常数，模拟 Rituximab 阻断效果

[t_exp, y_exp] = ode15s(@pathway_model_4_antibody_immune_cellular_humoral, ...
                        [0 500], x0_exp, options, para_exp);

%% ==========================
% 2. 绘制符合 Nature 要求的合并图
% ==========================

% 定义配色 (Colorblind-friendly / Nature Style)
% 蓝色: Control / Normal
c_norm = [0/255, 114/255, 189/255]; 
% 红色: Treatment / Rituximab (使用深红强调差异)
c_exp  = [162/255, 20/255, 47/255];

% 设置画布大小 (宽24cm, 高11cm) - 适当加宽以容纳大字体
figure('Units', 'centimeters', 'Position', [2, 2, 24, 11], 'Color', 'w');

% 使用 TiledLayout (1行2列)
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'normal');

% ===============================
% 子图 A: 细胞免疫 (Tc Cells)
% ===============================
ax1 = nexttile;
hold(ax1, 'on');

% 绘制曲线 (线条加粗到 2.5)
p1 = plot(ax1, t_norm, y_norm(:,4), 'Color', c_norm, 'LineWidth', 2.5);
p2 = plot(ax1, t_exp,  y_exp(:,4),  'Color', c_exp,  'LineWidth', 2.5);

% 设置坐标轴标签与标题
ylabel(ax1, 'Specific CD8^+ T Cells', 'FontWeight', 'bold', 'FontSize', 15);
xlabel(ax1, 'Time (AU)', 'FontWeight', 'bold', 'FontSize', 15);
title(ax1, 'Cellular Immunity', 'FontWeight', 'bold', 'FontSize', 16);

% 添加必要的标注 (Annotation)
% 在图中标注出实验组升高的现象，这是文章的核心论点
text(ax1, 250, max(y_exp(:,4))*0.9, {'\uparrow Enhanced Tc Response', '(Treatment Group)'}, ...
     'Color', c_exp, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% 设置通用属性
grid(ax1, 'on'); ax1.GridAlpha = 0.15;
set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Arial');

% 添加子图编号 'a' (固定在左上角)
text(ax1, -0.18, 1.08, 'a', 'Units', 'Normalized', 'FontSize', 22, 'FontWeight', 'bold', 'FontName', 'Arial');

% ===============================
% 子图 B: 体液免疫 (Antibody)
% ===============================
ax2 = nexttile;
hold(ax2, 'on');

% 绘制曲线 (Rituximab 组用虚线强调其受损/下降状态)
plot(ax2, t_norm, y_norm(:,1), 'Color', c_norm, 'LineWidth', 2.5);
plot(ax2, t_exp,  y_exp(:,1),  'Color', c_exp,  'LineWidth', 2.5, 'LineStyle', '--'); 

% 设置坐标轴标签与标题
ylabel(ax2, 'Specific Antibody Level', 'FontWeight', 'bold', 'FontSize', 15);
xlabel(ax2, 'Time (AU)', 'FontWeight', 'bold', 'FontSize', 15);
title(ax2, 'Humoral Immunity', 'FontWeight', 'bold', 'FontSize', 16);

% 添加标注，指出抗体受阻
text(ax2, 300, max(y_norm(:,1))*0.6, {'\downarrow Impaired Antibody', '(Treatment Group)'}, ...
     'Color', c_exp, 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% 设置通用属性
grid(ax2, 'on'); ax2.GridAlpha = 0.15;
set(ax2, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 14, 'FontName', 'Arial');

% 添加子图编号 'b'
text(ax2, -0.18, 1.08, 'b', 'Units', 'Normalized', 'FontSize', 22, 'FontWeight', 'bold', 'FontName', 'Arial');

% ===============================
% 全局图例 (Shared Legend)
% ===============================
% 在右侧子图中添加图例
lg = legend(ax2, {'Control Group (Unblocked)', 'Rituximab Group (Blocked Humoral)'}, ...
    'Location', 'northeast', 'FontSize', 12, 'EdgeColor', 'none'); 

hold off;