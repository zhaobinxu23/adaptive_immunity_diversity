function main_nature_antibody_dose()
    %% 1. 初始化设置
    clear; clc; close all;

    % 定义剂量等级数
    num_doses = 9;
    
    % --- Nature 风格配色: 灰红渐变 (Grey to Red) ---
    % 灰色代表低剂量/对照，红色越深代表剂量越高
    % 这种配色在肿瘤免疫学文章中非常经典
    start_color = [0.6, 0.6, 0.6]; % 灰色 (Control / Low)
    end_color   = [0.8, 0.0, 0.0]; % 深红 (High Dosage)
    
    % 生成颜色矩阵 (9行3列)
    colors = [linspace(start_color(1), end_color(1), num_doses)', ...
              linspace(start_color(2), end_color(2), num_doses)', ...
              linspace(start_color(3), end_color(3), num_doses)'];

    %% 2. 准备画布
    % 尺寸: 8.9 cm (Nature 单栏) x 7.5 cm
    figure('Units', 'centimeters', 'Position', [10, 10, 8.9, 7.5], 'Color', 'w');
    ax = axes;
    hold(ax, 'on');

    % 绘制 T=10 的垂直虚线 (表示给药时刻)
    xline(10, '--', 'Color', [0.3 0.3 0.3], 'LineWidth', 1, 'DisplayName', 'Treatment Start');

    %% 3. 循环仿真
    lines = gobjects(1, num_doses); % 存储句柄
    legend_str = cell(1, num_doses);

    % 基础参数定义 (保持原逻辑不变)
    % ---------------------------------
    base_x0 = zeros(7, 1);
    base_x0(4) = 1e5; % A
    base_x0(5) = 1e3; % Tc
    
    para = zeros(18, 1);
    para(1) = 1.1; para(2) = 0.05; para(3) = 10; para(4) = 1e7;
    para(5) = 1e-4; para(6) = 1e-7; para(7) = 0; para(8) = 1e4;
    para(9) = 1; para(10) = 2; para(11) = 1e-3; para(12) = 1e2;
    para(13) = 0.1; para(14) = 2.0e-4; para(15) = 1; para(16) = 1e-3;
    para(17) = 100; para(18) = 0.05;
    % ---------------------------------

    for i = 1:num_doses
        % --- 第一阶段: 0 -> 10 ---
        [t1, y1] = ode15s(@pathway_model_tumor_5_mono_antibody, [0 10], base_x0, [], para);
        
        % --- 计算第 10 秒的状态并添加抗体 ---
        current_x0 = zeros(7, 1);
        for k = 1:7
            current_x0(k) = interp1(t1, y1(:, k), 10);
        end
        
        % 注入抗体: Dosage 随 i 增加
        % i=1 -> 0, i=2 -> 1e8 ...
        dosage_val = (i-1) * 10^8;
        current_x0(6) = current_x0(6) + dosage_val;
        
        % --- 第二阶段: 10 -> 100 ---
        [t2, y2] = ode15s(@pathway_model_tumor_5_mono_antibody, [10 100], current_x0, [], para);
        
        % --- 拼接数据 (为了画出连续光滑的曲线) ---
        t_full = [t1; t2];
        y_full = [y1(:, 1); y2(:, 1)]; % 只取 tumor concentration (第1列)
        
        % --- 绘图 ---
        lines(i) = plot(t_full, y_full, ...
            'LineWidth', 1.5, ...       % 线宽 1.5 pt (比原代码 5 pt 更精致)
            'Color', colors(i, :));
        
        % 图例文字
        if i == 1
            legend_str{i} = 'Control (No Dose)';
        else
            % 使用科学计数法显示剂量
            legend_str{i} = sprintf('Dose: %d \\times 10^8', i-1);
        end
    end

    %% 4. 深度美化 (Nature Formatting)

    % --- 坐标轴 ---
    set(ax, 'Box', 'off', ...
            'TickDir', 'out', ...
            'LineWidth', 1, ...
            'XMinorTick', 'on', ...
            'YMinorTick', 'on', ...
            'FontName', 'Arial', ...
            'FontSize', 8, ...
            'XColor', [0.15 0.15 0.15], ...
            'YColor', [0.15 0.15 0.15]);

    % --- 标签 ---
    xlabel('Time (days)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel('Tumor Cell Concentration', 'FontName', 'Arial', 'FontSize', 9);

    % --- 范围与布局 ---
    xlim([0, 100]);
    % 自动调整 Y 轴上限，留白 10%
    all_y_data = cellfun(@(x) max(x), get(lines, 'YData')); 
    ylim([0, max(all_y_data) * 1.1]);

    % --- 图例优化 ---
    % 因为有 9 条线，全显示图例太占地方，我们只显示 "Control" 和 "Max Dose"
    % 或者把图例缩小放在旁边。这里演示精简图例法：
    subset_idx = [1, 3, 5, 7, 9]; % 只显示奇数编号的图例，保持整洁
    legend(lines(subset_idx), legend_str(subset_idx), ...
        'Location', 'best', ...
        'Box', 'off', ...
        'FontName', 'Arial', 'FontSize', 7);

    % 或者添加文字注释代替图例
    text(95, y_full(end), 'High Dose', ...
        'Color', end_color, 'FontSize', 7, 'FontName', 'Arial', ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');

    title('Tumor Response to Antibody Dosage', 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'normal');

    % print(gcf, 'Nature_Antibody_Response.tif', '-dtiff', '-r600');
    disp('绘图完成。');
end