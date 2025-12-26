%% 1. 数据检查与准备
if ~exist('y', 'var') || ~exist('y_new', 'var')
    error('请先运行仿真代码 generating t, y, t_new, y_new');
end

% 定义病毒索引 (根据你的描述)
v_idx = 2606; 

% 确定统一的 Y 轴范围 (至关重要，用于视觉对比)
% 获取两次模拟中的最大值和最小值
max_val = max(max(y(:, v_idx)), max(y_new(:, v_idx)));
min_val = 1e-2; % 设置一个合理的对数下限，避免 Log(0) 错误
common_ylim = [min_val, max_val * 5]; % 统一上限，留出标题空间

% 确定统一的 X 轴范围
common_xlim = [0, max(max(t), max(t_new))];

%% 2. 创建画布 (TiledLayout)
% 宽度 16cm (适合 Nature 双栏排版)，高度 8cm
figure('Units', 'centimeters', 'Position', [5, 5, 16, 8], 'Color', 'w');

% 创建 1x2 的分块布局 (左右并排)，紧凑间距
tlo = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% 3. Panel a: 初次感染 (Primary Infection / No Memory)
ax1 = nexttile;
% 使用半对数坐标 (Log Y)
p1 = semilogy(t, y(:, v_idx), 'LineWidth', 2.5);

% 样式设置
p1.Color = [0.85, 0.33, 0.1]; % 砖红色 (代表危险/严重)
title(ax1, 'Primary Infection (No Memory)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
ylabel(ax1, 'Viral Load (copies/mL)', 'FontName', 'Arial', 'FontSize', 9);
xlabel(ax1, 'Time (hours)', 'FontName', 'Arial', 'FontSize', 9);

% 关键：添加 'a' 标签
text(ax1, -0.15, 1.05, 'a', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

% 统一坐标轴
xlim(ax1, common_xlim);
ylim(ax1, common_ylim);

% 美化坐标轴
set(ax1, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 8, 'FontName', 'Arial');
grid(ax1, 'on'); ax1.GridAlpha = 0.1;

%% 4. Panel b: 二次感染 (Secondary Infection / With Memory)
ax2 = nexttile;
% 绘制曲线
p2 = semilogy(t_new, y_new(:, v_idx), 'LineWidth', 2.5);

% 样式设置
p2.Color = [0, 0.45, 0.74]; % 科学蓝 (代表控制/免疫)
title(ax2, 'Re-infection (With Memory)', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');

% 为了排版整洁，右图可以省略 Y 轴标签 (共用左边的)，或者保留
% 这里我们保留但留空，或者设为 YTickLabel 为空以节省空间
% set(ax2, 'YTickLabel', []); 
xlabel(ax2, 'Time (hours)', 'FontName', 'Arial', 'FontSize', 9);

% 关键：添加 'b' 标签
text(ax2, -0.15, 1.05, 'b', 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 12, 'FontWeight', 'bold');

% 在图中标注具体的抑制倍数 (Nature 喜欢的定量描述)
peak_1 = max(y(:, v_idx));
peak_2 = max(y_new(:, v_idx));
fold_change = peak_1 / peak_2;
dim_str = sprintf('Fold reduction: ~10^{%.1f}', log10(fold_change));
% 在图内添加文本框
text(ax2, 0.05, 0.9, dim_str, 'Units', 'normalized', 'FontName', 'Arial', 'FontSize', 8, 'EdgeColor', 'none', 'BackgroundColor', 'none');

% 统一坐标轴
xlim(ax2, common_xlim);
ylim(ax2, common_ylim);

% 美化坐标轴
set(ax2, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 8, 'FontName', 'Arial');
grid(ax2, 'on'); ax2.GridAlpha = 0.1;

%% 5. 导出建议
% 导出为高分辨率 TIFF
% print(gcf, 'Figure_Comparison_Merged.tif', '-dtiff', '-r600');