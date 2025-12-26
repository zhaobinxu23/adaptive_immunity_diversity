function main_nature_neoantigen_timing()
    %% 1. 初始化设置
    clear; clc; close all;

    % 定义两个给药时间点
    injection_times = [5, 10]; 
    panel_titles = {'Early Intervention (Day 5)', 'Late Intervention (Day 10)'};
    
    % 剂量设置 (0 到 8 * 10^8)
    num_doses = 9;
    
    % --- Nature 配色: 渐变红 (Sequential Red Palette) ---
    % 颜色越深，剂量越高
    c_start = [1.0, 0.85, 0.80]; % 极浅红 (Control/Low)
    c_end   = [0.7, 0.0, 0.0];   % 深红 (High Dosage)
    colors = [linspace(c_start(1), c_end(1), num_doses)', ...
              linspace(c_start(2), c_end(2), num_doses)', ...
              linspace(c_start(3), c_end(3), num_doses)'];

    %% 2. 准备画布 (双栏布局)
    % 尺寸: 18 cm (Nature 双栏宽度) x 8 cm
    figure('Units', 'centimeters', 'Position', [5, 5, 18, 8], 'Color', 'w');

    % 基础参数 (保持你原代码的参数)
    base_x0 = [0; 0; 0; 1e5; 1e3]; % Tu, V, C, A, Tc
    para = zeros(18, 1);
    para(1)=1.1; para(2)=0.05; para(3)=10; para(4)=1e7; para(5)=1e-4;
    para(6)=1e-7; para(7)=0; para(8)=1e4; para(9)=1; para(10)=2;
    para(11)=1e-3; para(12)=1e2; para(13)=0.1; para(14)=2.0e-4;
    para(15)=1; para(16)=1e-3; para(17)=100; para(18)=0.05;

    %% 3. 主循环：分别绘制左图 (T=5) 和右图 (T=10)
    for p = 1:2
        % 创建子图
        ax = subplot(1, 2, p);
        hold(ax, 'on');
        
        T_inject = injection_times(p);
        lines = gobjects(1, num_doses);
        
        % 绘制给药时间点的竖虚线
        xline(T_inject, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, ...
              'HandleVisibility', 'off'); % 不在图例显示
        
        % 文字标注给药时刻
        text(T_inject, 0, sprintf(' Tx (Day %d)', T_inject), ...
             'VerticalAlignment', 'bottom', 'FontSize', 7, 'FontName', 'Arial', ...
             'Color', [0.3 0.3 0.3]);

        % --- 剂量循环 ---
        for i = 1:num_doses
            dosage = (i-1) * 10^8;
            
            % Phase 1: 0 -> T_inject
            [t1, y1] = ode15s(@pathway_model_tumor_5, [0 T_inject], base_x0, [], para);
            
            % 计算中间状态 (Interp1)
            curr_x0 = zeros(5, 1);
            for k = 1:5
                curr_x0(k) = interp1(t1, y1(:, k), T_inject);
            end
            
            % *** 施加干预 *** (对应原代码: x0(2) + dosage)
            curr_x0(2) = curr_x0(2) + dosage;
            
            % Phase 2: T_inject -> 100
            [t2, y2] = ode15s(@pathway_model_tumor_5, [T_inject 100], curr_x0, [], para);
            
            % 数据拼接 (Stitch Data)
            t_full = [t1; t2];
            y_full = [y1(:, 1); y2(:, 1)]; % 仅提取 Tumor Concentration
            
            % 绘图
            lines(i) = plot(t_full, y_full, ...
                'LineWidth', 2, ...         % 2pt 线宽，清晰且不臃肿
                'Color', colors(i, :));
        end
        
        %% 4. 子图美化 (Formatting)
        set(ax, 'Box', 'off', ...
                'TickDir', 'out', ...
                'LineWidth', 1, ...
                'XMinorTick', 'on', ...
                'YMinorTick', 'on', ...
                'FontName', 'Arial', ...
                'FontSize', 8, ...
                'XColor', [0.1 0.1 0.1], ...
                'YColor', [0.1 0.1 0.1]);
        
        % 统一坐标轴范围 (便于左右对比)
        xlim([0, 100]);
        % 根据数据自动调整Y轴，这里假设最大值统一
        ylim([0, max(y_full)*1.1]); 

        title(panel_titles{p}, 'FontName', 'Arial', 'FontSize', 10, 'FontWeight', 'bold');
        xlabel('Time (days)', 'FontName', 'Arial', 'FontSize', 9);
        
        if p == 1
            ylabel('Tumor Concentration', 'FontName', 'Arial', 'FontSize', 9);
            % 左上角添加 "a" 标签 (Nature 习惯)
            text(-0.15, 1.05, 'a', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
        else
             % 右上角添加 "b" 标签
            text(-0.15, 1.05, 'b', 'Units', 'normalized', 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial');
            
            % 仅在右图添加图例 (精简显示)
            % 只显示 Control, 中间量, 最大量
            idx_show = [1, 5, 9];
            labels = {'Control', 'Med Dose', 'High Dose'};
            legend(lines(idx_show), labels, ...
                   'Location', 'northeast', 'Box', 'off', 'FontSize', 7, 'FontName', 'Arial');
        end
    end

    % 导出建议
    % print(gcf, 'Fig_Neoantigen_Timing_Comparison.tif', '-dtiff', '-r600');
    disp('Nature 风格绘图完成。');
end