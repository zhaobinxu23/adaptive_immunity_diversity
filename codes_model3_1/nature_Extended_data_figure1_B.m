figure('Color', 'w', 'Position', [100, 100, 1400, 800]);
    hold on;
    axis equal;
    
    % ==========================================================
    % 1. 数据准备
    % ==========================================================
    
    % 原始相对数值
    raw_str = [...
        '8.21691236608131e-14 2.52627846604407e-11 1.75466537269336e-09 2.85046633638077e-08 1.13041112249312e-07 1.13041112249312e-07 2.85046633638077e-08 1.75466537269336e-09 2.52627846604407e-11 8.21691236608131e-14 ' ...
        '2.52627846604407e-11 7.76700858383578e-09 5.39469468415719e-07 8.76371975622824e-06 3.47543353184500e-05 3.47543353184500e-05 8.76371975622824e-06 5.39469468415719e-07 7.76700858383578e-09 2.52627846604407e-11 ' ...
        '1.75466537269336e-09 5.39469468415719e-07 3.74696775742474e-05 0.000608697568337430 0.00241391555024220 0.00241391555024220 0.000608697568337430 3.74696775742474e-05 5.39469468415719e-07 1.75466537269336e-09 ' ...
        '2.85046633638077e-08 8.76371975622824e-06 0.000608697568337430 0.00988833514688556 0.0392142292308970 0.0392142292308970 0.00988833514688556 0.000608697568337430 8.76371975622824e-06 2.85046633638077e-08 ' ...
        '1.13041112249312e-07 3.47543353184500e-05 0.00241391555024220 0.0392142292308970 0.155512101009002 0.155512101009002 0.0392142292308970 0.00241391555024220 3.47543353184500e-05 1.13041112249312e-07 ' ...
        '1.13041112249312e-07 3.47543353184500e-05 0.00241391555024220 0.0392142292308970 0.155512101009002 0.155512101009002 0.0392142292308970 0.00241391555024220 3.47543353184500e-05 1.13041112249312e-07 ' ...
        '2.85046633638077e-08 8.76371975622824e-06 0.000608697568337430 0.00988833514688556 0.0392142292308970 0.0392142292308970 0.00988833514688556 0.000608697568337430 8.76371975622824e-06 2.85046633638077e-08 ' ...
        '1.75466537269336e-09 5.39469468415719e-07 3.74696775742474e-05 0.000608697568337430 0.00241391555024220 0.00241391555024220 0.000608697568337430 3.74696775742474e-05 5.39469468415719e-07 1.75466537269336e-09 ' ...
        '2.52627846604407e-11 7.76700858383578e-09 5.39469468415719e-07 8.76371975622824e-06 3.47543353184500e-05 3.47543353184500e-05 8.76371975622824e-06 5.39469468415719e-07 7.76700858383578e-09 2.52627846604407e-11 ' ...
        '8.21691236608131e-14 2.52627846604407e-11 1.75466537269336e-09 2.85046633638077e-08 1.13041112249312e-07 1.13041112249312e-07 2.85046633638077e-08 1.75466537269336e-09 2.52627846604407e-11 8.21691236608131e-14'];
    
    AA_raw = str2num(raw_str); 
    
    % 绝对数量计算
    TOTAL_POOL = 4 * 10^19;
    AA = AA_raw * TOTAL_POOL;
    
    % ==========================================================
    % 2. 绘图数据构建
    % ==========================================================
    plot_data = [];
    for m = 1:10
        for n = 1:10
            idx = 10*(m-1) + n;
            val = AA(idx);
            x_val = -(n - m) - 22; 
            plot_data = [plot_data; idx, x_val, val];
        end
    end
    
    % 排序保证堆叠
    [~, sort_idx] = sortrows(plot_data, [2, 1]);
    plot_data = plot_data(sort_idx, :);
    
    % ==========================================================
    % 3. 颜色映射 (Jet)
    % ==========================================================
    colormap(jet(256));
    c_map = colormap;
    
    min_val = min(plot_data(:,3));
    max_val = max(plot_data(:,3));
    min_log = log10(min_val);
    max_log = log10(max_val);
    
    % ==========================================================
    % 4. 绘图执行
    % ==========================================================
    radius = 0.35;
    
    count_map = containers.Map('KeyType','double','ValueType','double');
    total_map = containers.Map('KeyType','double','ValueType','double');
    
    % 统计每列总数
    for i = 1:size(plot_data, 1)
        x = plot_data(i, 2);
        if isKey(total_map, x)
            total_map(x) = total_map(x) + 1;
        else
            total_map(x) = 1;
            count_map(x) = 0;
        end
    end
    
    for i = 1:size(plot_data, 1)
        id = plot_data(i, 1);
        x = plot_data(i, 2);
        val = plot_data(i, 3);
        
        % Y坐标
        count_map(x) = count_map(x) + 1;
        y = count_map(x) - 1 - (total_map(x) - 1) / 2.0;
        
        % 颜色计算
        log_val = log10(val);
        color_idx = floor( (log_val - min_log) / (max_log - min_log) * 255 ) + 1;
        if color_idx > 256, color_idx = 256; end
        if color_idx < 1, color_idx = 1; end
        this_color = c_map(color_idx, :);
        
        % 绘制圆形
        rectangle('Position', [x-radius, y-radius, 2*radius, 2*radius], ...
                  'Curvature', [1 1], ...
                  'FaceColor', this_color, ...
                  'EdgeColor', [0.3 0.3 0.3], 'LineWidth', 0.5);
              
        % 文字颜色逻辑
        ratio = (color_idx - 1) / 255;
        if ratio < 0.35 || ratio > 0.85
            txt_color = 'w';
        else
            txt_color = 'k';
        end
        
        label_str = sprintf('Ab%d', id);
        
        % 【关键修改】FontSize 增大至 7
        text(x, y, label_str, ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'FontSize', 7, ...  % 这里调大字号
             'FontWeight', 'bold', ...
             'Color', txt_color);
    end
    
    % ==========================================================
    % 5. 坐标轴与Colorbar
    % ==========================================================
    x_unique = unique(plot_data(:,2));
    xlim([-32.5, -11.5]);
    ylim([-6, 6]);
    xticks(min(x_unique):1:max(x_unique));
    xlabel('Affinity Group (log_{10} K_d)', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'YTick', [], 'YColor', 'none', 'Box', 'off');
    
    text(-31, -6.5, '\leftarrow Weak Binding', 'Color', [0.6 0 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
    text(-13, -6.5, 'Strong Binding \rightarrow', 'Color', [0 0.5 0], 'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment','center');
    
    title('Initial Antibody Absolute Number Distribution (N_0)', 'FontSize', 14, 'FontWeight', 'bold');
    
    c = colorbar;
    c.Label.String = 'Absolute Number (N_0)';
    c.Label.FontSize = 11;
    c.Label.FontWeight = 'bold';
    caxis([min_log, max_log]);
    
    ticks_log = ceil(min_log):2:floor(max_log); 
    c.Ticks = ticks_log;
    
    labels = cell(length(ticks_log), 1);
    for i = 1:length(ticks_log)
        labels{i} = sprintf('10^{%d}', ticks_log(i));
    end
    c.TickLabels = labels;
