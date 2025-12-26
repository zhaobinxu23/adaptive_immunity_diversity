function main_nature_figure5_g()
    %% 1. 初始化与设置
    clear; clc; close all;

    % --- 实验参数设置 ---
    % 定义要扫描的参数 para(10) 的值
    % 原代码变化: 2.0 -> 2.1 -> 2.2
    param_values = [2.0, 2.1, 2.2];
    num_trials = length(param_values);
    
    % --- Nature 风格配色 (NPG Style) ---
    % 使用红/蓝/绿三种高对比度颜色，分别代表参数的变化
    colors = [
        0.25, 0.60, 0.81; % 蓝色 (2.0 - 基准)
        0.00, 0.62, 0.45; % 绿色 (2.1 - 中等)
        0.89, 0.10, 0.11  % 红色 (2.2 - 高)
    ];

    %% 2. 准备画布
    % 尺寸: 8.9 cm (Nature 单栏标准宽) x 7 cm
    figure('Units', 'centimeters', 'Position', [10, 10, 8.9, 7], 'Color', 'w');
    ax = axes;
    hold(ax, 'on');

    %% 3. 模型基础参数定义
    % 初始条件
    x0 = zeros(5, 1);
    x0(1)=0;    % Tu
    x0(2)=0;    % V
    x0(3)=0;    % C
    x0(4)=1e5;  % A
    x0(5)=1e3;  % Tc

    % 此处定义基础参数
    para = zeros(18, 1);
    para(1)=1.1; para(2)=0.05; para(3)=10; para(4)=1e7;
    para(5)=1e-4; para(6)=1e-7; para(7)=0; para(8)=1e4;
    para(9)=1; 
    % para(10) 将在循环中赋值
    para(11)=1e-3; para(12)=1e2; para(13)=0.1; para(14)=2.0e-4;
    para(15)=1; para(16)=1e-3; para(17)=100; para(18)=0.05;

    %% 4. 循环仿真
    lines = gobjects(1, num_trials);
    legend_str = cell(1, num_trials);

    for i = 1:num_trials
        % 更新参数 para(10)
        para(10) = param_values(i);
        
        % 运行 ODE
        [t, y] = ode15s(@pathway_model_tumor_5, [0 100], x0, [], para);
        
        % 绘图
        lines(i) = plot(t, y(:, 1), ...
            'LineWidth', 2, ...          % 线宽改为 2 (精致且清晰)
            'Color', colors(i, :));      % 使用 NPG 配色
        
        % 记录图例 (使用 LaTeX 格式)
        legend_str{i} = sprintf('k_4 = %.1f', param_values(i));
    end

    %% 5. 深度美化 (Nature Formatting)

    % --- 坐标轴样式 ---
    set(ax, 'Box', 'off', ...           % 去掉右侧和顶部的框线
            'TickDir', 'out', ...       % 刻度向外
            'LineWidth', 1, ...         % 坐标轴线宽 1pt
            'XMinorTick', 'on', ...     % 开启次级刻度
            'YMinorTick', 'on', ...
            'FontName', 'Arial', ...    % 强制使用 Arial 字体
            'FontSize', 8, ...          % 字号 8pt
            'XColor', [0.15 0.15 0.15], ... 
            'YColor', [0.15 0.15 0.15]);

    % --- 标签与标题 ---
    xlabel('Time (arbitrary units)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel('Tumor Concentration', 'FontName', 'Arial', 'FontSize', 9);

    % --- 范围微调 ---
    xlim([0, 100]);
    % 自动获取最大 Y 值并留出 10% 顶空，防止曲线顶到头
    ymax = max(cellfun(@max, get(lines, 'YData')));
    ylim([0, ymax * 1.1]);

    % --- 图例优化 ---
    legend(lines, legend_str, ...
        'Location', 'northwest', ...    % 根据曲线位置调整，通常放左上
        'Box', 'off', ...               % 去掉图例边框
        'FontName', 'Arial', 'FontSize', 7);

    % --- 导出建议 ---
    % print(gcf, 'Nature_Parameter_Study.tif', '-dtiff', '-r600');
    disp('绘图完成。');
end