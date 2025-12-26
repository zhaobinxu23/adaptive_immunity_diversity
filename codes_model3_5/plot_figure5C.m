function main_pd1_nature_plot()
    %% 1. 初始化设置
    clear; clc; close all;
    
    % --- 定义 4 种实验条件 (例如：PD-1 抗体剂量/效率) ---
    % 0 = 对照组(无治疗), 逐步增加治疗强度
    conditions = [0, 0.3, 0.6, 0.9]; 
    cond_names = {'Control (No Tx)', 'Low Dose', 'Medium Dose', 'High Dose'};
    
    % 时间跨度
    t_span = [0, 60]; % 天
    y0 = [100; 10; 10]; % 初始值: [肿瘤, 效应T, 耗竭T] (根据你的模型调整)

    %% 2. 准备画布 (Nature Standard)
    % 单栏插图标准宽度: 8.5 - 9.0 cm
    figure('Units', 'centimeters', 'Position', [10, 10, 9, 7.5], 'Color', 'w');
    ax = axes;
    hold(ax, 'on');

    %% 3. 定义 Nature 专属配色 (NPG Palette)
    % 这种配色对比度高，能够区分 Control 和不同程度的治疗组
    % 顺序：砖红(对照), 蓝色, 绿色, 紫色
    colors = [
        0.8500, 0.3250, 0.0980; % Vermilion (Control)
        0.0000, 0.4470, 0.7410; % Blue
        0.4660, 0.6740, 0.1880; % Green
        0.4940, 0.1840, 0.5560  % Purple
    ];

    %% 4. 循环仿真与绘图
    lines = gobjects(1, 4); % 存储句柄用于图例
    
    for i = 1:4
        % 获取当前条件参数
        curr_effect = conditions(i);
        
        % --- 运行 ODE 求解 ---
        % 注意：这里调用底部的 pd1_ode 函数
        [t, y] = ode45(@(t, y) pd1_ode(t, y, curr_effect), t_span, y0);
        
        % --- 绘制 y(:,1) 肿瘤曲线 ---
        lines(i) = plot(t, y(:, 1), ...
            'LineWidth', 2, ...         % 线宽 2磅，清晰
            'Color', colors(i, :));     % 分配颜色
    end

    %% 5. 美化细节 (Nature Style Formatting)
    
    % --- 坐标轴 ---
    set(ax, 'Box', 'off', ...           % 去掉右/上边框
            'TickDir', 'out', ...       % 刻度朝外
            'LineWidth', 1, ...         % 坐标轴线宽合适
            'XMinorTick', 'on', ...     % 开启次级刻度
            'YMinorTick', 'on', ...
            'FontName', 'Arial', ...    % 必须用 Arial/Helvetica
            'FontSize', 8, ...          % 8pt 标准字号
            'XColor', [0.1, 0.1, 0.1], ... % 深灰色坐标轴
            'YColor', [0.1, 0.1, 0.1]);

    % --- 标签 ---
    xlabel('Time (days)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel('Tumor Volume (mm^3)', 'FontName', 'Arial', 'FontSize', 9);
    
    % --- 范围微调 ---
    xlim([0, 60]);
    ylim([0, max(ylim)*1.05]); % 稍微留出顶部空间

    % --- 图例 ---
    % 无边框，字体略小
    legend(lines, cond_names, ...
        'Location', 'best', ...
        'Box', 'off', ...
        'FontName', 'Arial', 'FontSize', 7);

    % --- 导出建议 ---
    % print(gcf, 'PD1_Effect_Comparison.tif', '-dtiff', '-r600');
    
    disp('绘图完成。');
end

%% --- 附件中的核心模型方程 (需根据你实际方程修改) ---
function dydt = pd1_ode(t, y, effect_param)
    % 假设 y(1) = 肿瘤, y(2) = T细胞
    % effect_param = PD-1 阻断带来的 T 细胞杀伤增强系数
    
    T = y(1); % Tumor
    E = y(2); % Effector Cells
    
    % 简单的 Logistic 生长 + 免疫杀伤模型
    r = 0.5;      % 肿瘤生长率
    K = 3000;     % 环境容纳量
    kill_base = 0.05; % 基础杀伤率
    
    % 加入参数影响：efficacy 越高，杀伤越强
    killing_rate = kill_base * (1 + effect_param * 5); 
    
    dT = r * T * (1 - T/K) - killing_rate * E * T;
    dE = 0.1 * T - 0.2 * E; % 简单的 T 细胞募集与衰减
    
    % 这里的维度必须和你 ode45 初始值对应
    dydt = [dT; dE; 0]; 
end