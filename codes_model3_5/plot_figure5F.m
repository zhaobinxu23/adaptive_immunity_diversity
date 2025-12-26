function main_nature_figure5_f()
    %% 1. 初始化设置
    clear; clc; close all;

    % --- 实验设置 ---
    % 定义要比较的 k5 参数值 (从低到高)
    k5_values = [2.0e-4, 4.0e-4, 6.0e-4];
    num_trials = length(k5_values);
    
    % --- Nature 风格排版 ---
    % 画布大小：8.9cm (单栏宽度) x 7cm
    figure('Units', 'centimeters', 'Position', [10, 10, 8.9, 7], 'Color', 'w');
    ax = axes;
    hold(ax, 'on');

    % --- NPG 风格配色 (Blue, Green, Red) ---
    % 这种配色对比强烈，且符合色盲友好原则
    colors = [
        0.25, 0.60, 0.81; % #4DBBD5 (Blue - Low k5)
        0.00, 0.62, 0.45; % #00A087 (Green - Medium k5)
        0.89, 0.10, 0.11  % #E64B35 (Red - High k5)
    ];

    %% 2. 模型参数定义 (保持原逻辑)
    % 初始条件
    x0 = zeros(5, 1);
    x0(1)=0; x0(2)=0; x0(3)=0; x0(4)=1e5; x0(5)=1e3;
    
    % 基础参数
    para = zeros(18, 1);
    para(1)=1.1; para(2)=0.05; para(3)=10; para(4)=1e7; para(5)=1e-4;
    para(6)=1e-7; para(7)=0; para(8)=1e4; para(9)=1; para(10)=2;
    para(11)=1e-3; para(12)=1e2; para(13)=0.1;
    % para(14) 会在循环中改变
    para(15)=1; para(16)=1e-3; para(17)=100; para(18)=0.05;

    %% 3. 循环仿真与绘图
    lines = gobjects(1, num_trials); % 存储句柄
    legend_str = cell(1, num_trials);
    
    for i = 1:num_trials
        % 更新 k5 参数
        para(14) = k5_values(i);
        
        % 运行求解器
        [t, y] = ode15s(@pathway_model_tumor_5, [0 100], x0, [], para);
        
        % 绘图
        lines(i) = plot(t, y(:,1), ...
            'LineWidth', 2, ...         % 标准线宽 1.5 - 2 pt
            'Color', colors(i, :));
        
        % 生成图例文字 (使用 LaTeX 格式)
        % 将 2.0e-4 格式化为 2.0 x 10^{-4}
        val_exp = log10(k5_values(i));
        val_base = k5_values(i) / 10^floor(val_exp);
        legend_str{i} = sprintf('k_5 = %.1f \\times 10^{%d}', val_base, floor(val_exp));
    end

    %% 4. 深度美化 (Nature Formatting)
    
    % --- 坐标轴设置 ---
    set(ax, 'Box', 'off', ...           % 去除右/上边框
            'TickDir', 'out', ...       % 刻度朝外
            'LineWidth', 1, ...         % 坐标轴线宽 1pt
            'XMinorTick', 'on', ...     % 开启次刻度
            'YMinorTick', 'on', ...
            'FontName', 'Arial', ...    % 强制使用 Arial
            'FontSize', 8, ...          % 字体大小 8pt
            'XColor', [0.15 0.15 0.15], ... % 深灰而非纯黑，视觉更柔和
            'YColor', [0.15 0.15 0.15]);

    % --- 标签 ---
    xlabel('Time (days)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel('Tumor Cell Concentration', 'FontName', 'Arial', 'FontSize', 9);
    
    % --- 范围微调 ---
    xlim([0, 100]);
    % 动态获取最大Y值并留出 10% 顶空
    ymax = max(cellfun(@max, get(lines, 'YData')));
    ylim([0, ymax * 1.1]);

    % --- 图例 ---
    % 显示在右上角，无边框
    legend(lines, legend_str, ...
        'Location', 'northeast', ...
        'Box', 'off', ...
        'FontName', 'Arial', 'FontSize', 7, ...
        'Interpreter', 'tex'); % 启用 TeX 解析数学符号

    % --- 标题 (可选，Nature 图表通常不需要标题，而是写在 Figure Legend 中) ---
    % 如果需要保留：
    % title('Effect of Parameter k_5 on Tumor Growth', 'FontSize', 10, 'FontWeight', 'normal');

    % print(gcf, 'Fig_k5_Comparison.tif', '-dtiff', '-r600');
    disp('Nature 风格绘图完成。');
end