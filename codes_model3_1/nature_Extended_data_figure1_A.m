 figure('Color', 'w', 'Position', [100, 100, 1400, 800]);
    ax_main = axes('Position', [0.05, 0.1, 0.78, 0.8]);
    hold(ax_main, 'on');
    axis(ax_main, 'equal');

    % ==========================================================
    % 1. 参数定义与数据准备
    % ==========================================================
    
    % k_on 范围: 5*10^-22 到 5*10^-13
    min_log_kon = log10(5 * 10^(-22)); 
    max_log_kon = log10(5 * 10^(-13)); 
    
    % k_off 范围: 10^0 到 10^9
    min_log_koff = 0; 
    max_log_koff = 9; 
    
    % 准备自定义色谱 (Blues 和 Reds)
    blues_map = [linspace(0.95, 0.05, 256)', linspace(0.95, 0.25, 256)', linspace(1, 0.6, 256)'];
    reds_map  = [linspace(1, 0.4, 256)', linspace(0.95, 0, 256)', linspace(0.95, 0, 256)'];

    data_struct = [];
    
    for m = 1:10
        for n = 1:10
            ab_id = 10 * (m - 1) + n;
            val_kon = 5 * 10^(m - 23);
            val_koff = 10^(n - 1);
            log_kon_val = log10(val_kon);
            log_koff_val = log10(val_koff);
            x_val = -(n - m) - 22;
            data_struct = [data_struct; ab_id, x_val, log_kon_val, log_koff_val, m];
        end
    end
    
    [~, sort_idx] = sortrows(data_struct, [2, 5]);
    data_struct = data_struct(sort_idx, :);

    % ==========================================================
    % 2. 绘图循环 (绘制半圆)
    % ==========================================================
    
    radius = 0.35;
    x_counts = containers.Map('KeyType','double','ValueType','double');
    x_totals = containers.Map('KeyType','double','ValueType','double');
    
    for i = 1:size(data_struct, 1)
        x = data_struct(i, 2);
        if isKey(x_totals, x)
            x_totals(x) = x_totals(x) + 1;
        else
            x_totals(x) = 1;
            x_counts(x) = 0;
        end
    end
    
    theta_left = linspace(pi/2, 3*pi/2, 50);
    theta_right = linspace(-pi/2, pi/2, 50);
    
    for i = 1:size(data_struct, 1)
        id = data_struct(i, 1);
        x = data_struct(i, 2);
        l_kon = data_struct(i, 3);
        l_koff = data_struct(i, 4);
        
        x_counts(x) = x_counts(x) + 1;
        y = x_counts(x) - 1 - (x_totals(x) - 1) / 2.0;
        
        % --- 左半圆 (Kon - Blue) ---
        norm_idx = floor( (l_kon - min_log_kon) / (max_log_kon - min_log_kon) * 255 ) + 1;
        norm_idx = max(1, min(256, norm_idx));
        c_left = blues_map(norm_idx, :);
        
        X_L = x + radius * cos(theta_left);
        Y_L = y + radius * sin(theta_left);
        
        patch(ax_main, 'XData', X_L, 'YData', Y_L, 'FaceColor', c_left, 'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.1);
        
        % --- 右半圆 (Koff - Red) ---
        norm_idx_off = floor( (l_koff - min_log_koff) / (max_log_koff - min_log_koff) * 255 ) + 1;
        norm_idx_off = max(1, min(256, norm_idx_off));
        c_right = reds_map(norm_idx_off, :);
        
        X_R = x + radius * cos(theta_right);
        Y_R = y + radius * sin(theta_right);
        
        patch(ax_main, 'XData', X_R, 'YData', Y_R, 'FaceColor', c_right, 'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.1);
        
        % 文字标签设置
        if norm_idx > 150 || norm_idx_off > 150 
             txt_col = 'w';
        else
             txt_col = 'k';
        end
        if strcmp(txt_col, 'k'), txt_col = [0.1 0.1 0.1]; end
        
        % 【修改点】FontSize 调整为 7
        text(ax_main, x, y, sprintf('Ab\n%d', id), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'FontSize', 7, ...          % 增大字体
             'FontWeight', 'bold', ...
             'Color', txt_col);
    end
    
    % ==========================================================
    % 3. 主坐标轴修饰
    % ==========================================================
    unique_x = sort(cell2mat(keys(x_totals)));
    xlim(ax_main, [-32.5, -11.5]);
    ylim(ax_main, [-6, 6]);
    xticks(ax_main, unique_x);
    xlabel(ax_main, 'Affinity Group (log_{10} K_d)', 'FontSize', 12, 'FontWeight', 'bold');
    
    set(ax_main, 'YTick', [], 'YColor', 'none', 'Box', 'off');
    
    text(ax_main, -31, -6.5, '\leftarrow Weak Binding', 'Color', [0.6 0 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    text(ax_main, -13, -6.5, 'Strong Binding \rightarrow', 'Color', [0 0.5 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    
    title(ax_main, 'IgM Affinity Clusters: 100 Antibodies Mapped by Kinetic Rates', 'FontSize', 14, 'FontWeight', 'bold');

    % ==========================================================
    % 4. 绘制双 Colorbar
    % ==========================================================
    ax_cb1 = axes('Position', [0.85, 0.15, 0.025, 0.3], 'Visible', 'off');
    colormap(ax_cb1, blues_map);
    caxis(ax_cb1, [min_log_kon max_log_kon]); 
    cb1 = colorbar(ax_cb1, 'Position', [0.85, 0.15, 0.025, 0.3]);
    cb1.Label.String = 'k_{on} Value (M^{-1} s^{-1})';
    cb1.Label.FontSize = 10;
    cb1.Label.FontWeight = 'bold';
    
    ticks_val_kon = ceil(min_log_kon):2:floor(max_log_kon);
    cb1.Ticks = ticks_val_kon;
    labels_kon = cell(length(ticks_val_kon), 1);
    for k = 1:length(ticks_val_kon)
        labels_kon{k} = sprintf('10^{%d}', round(ticks_val_kon(k))); 
    end
    cb1.TickLabels = labels_kon;
    
    ax_cb2 = axes('Position', [0.85, 0.55, 0.025, 0.3], 'Visible', 'off');
    colormap(ax_cb2, reds_map);
    caxis(ax_cb2, [min_log_koff max_log_koff]);
    cb2 = colorbar(ax_cb2, 'Position', [0.85, 0.55, 0.025, 0.3]);
    cb2.Label.String = 'k_{off} Value (s^{-1})';
    cb2.Label.FontSize = 10;
    cb2.Label.FontWeight = 'bold';
    
    ticks_val_koff = 0:1:9;
    cb2.Ticks = ticks_val_koff;
    labels_koff = cell(length(ticks_val_koff), 1);
    for k = 1:length(ticks_val_koff)
        labels_koff{k} = sprintf('10^{%d}', ticks_val_koff(k));
    end
    cb2.TickLabels = labels_koff;
    
    % ==========================================================
    % 5. 图例说明结构 (右上角小图)
    % ==========================================================
    ax_legend = axes('Position', [0.82, 0.88, 0.15, 0.1], 'Visible', 'off');
    hold(ax_legend, 'on');
    xlim(ax_legend, [0 1]); ylim(ax_legend, [0 1]);
    
    demo_r = 0.2;
    dx = 0.3; dy = 0.5;
    
    % 使用修复后的 patch 语法
    patch('Parent', ax_legend, ...
          'XData', dx + demo_r*cos(theta_left), ...
          'YData', dy + demo_r*sin(theta_left), ...
          'FaceColor', blues_map(180,:), ...
          'EdgeColor', 'gray');
          
    patch('Parent', ax_legend, ...
          'XData', dx + demo_r*cos(theta_right), ...
          'YData', dy + demo_r*sin(theta_right), ...
          'FaceColor', reds_map(180,:), ...
          'EdgeColor', 'gray');
        
    text(ax_legend, 0.6, 0.6, 'Node Structure', 'FontWeight', 'bold', 'FontSize', 10);
    text(ax_legend, 0.3, 0.2, 'left: k_{on} | right: k_{off}', 'HorizontalAlignment','center', 'FontSize', 8);