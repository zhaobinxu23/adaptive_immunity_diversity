function nature_figure1I_fusion()
    
    % 设置全局绘图属性
    set(0, 'DefaultAxesFontName', 'Arial');
    set(0, 'DefaultAxesFontSize', 10);
    set(0, 'DefaultLineLineWidth', 1.5);

    %% 1. 数据准备
    % 检查是否存在原始ODE文件，否则生成模拟数据
    if exist('pathway_model_many_antibody_immune_include_plasma', 'file') == 2
        fprintf('检测到模型文件，正在运行计算...\n');
        [time_points, values_kd, freq_IgM, freq_IgG] = run_simulation();
    else
        warning('未检测到ODE文件，正在生成演示用模拟数据...');
        [time_points, values_kd, freq_IgM, freq_IgG] = generate_mock_data();
    end

    %% 2. 绘图布局初始化
    figure('Units', 'pixels', 'Position', [50, 50, 1200, 900], 'Color', 'w');
    % 使用 3行 x 6列 的布局以实现灵活对齐
    t = tiledlayout(3, 6, 'TileSpacing', 'compact', 'Padding', 'compact');

    % 颜色方案：灰(T=0), 红(T=50, 爆发), 蓝(T=1000, 记忆)
    target_times = [0, 50, 1000];
    colors = {[0.5 0.5 0.5], [0.85 0.33 0.1], [0 0.45 0.74]}; 

    %% 3. 第一行：三维全景图 (3D Landscape)
    % 左侧：IgM (占据第一行前3列)
    ax_3d_m = nexttile(1, [1 3]);
    plot_3d_surface(ax_3d_m, time_points, values_kd, freq_IgM, target_times, colors, 'IgM');

    % 右侧：IgG (占据第一行后3列)
    ax_3d_g = nexttile(4, [1 3]);
    plot_3d_surface(ax_3d_g, time_points, values_kd, freq_IgG, target_times, colors, 'IgG');

    %% 4. 第二行 & 第三行：流式细胞图 (Flow Cytometry Snapshots)
    
    % --- Row 2: IgM 流式图 ---
    for i = 1:3
        % 计算位置：占据第2行的位置
        tile_idx = 6 + (i-1)*2 + 1; % 7, 9, 11
        ax = nexttile(tile_idx, [1 2]);
        
        t_val = target_times(i);
        c_code = colors{i};
        
        % 提取对应时刻数据
        t_idx = find(time_points == t_val, 1);
        if isempty(t_idx), t_idx=1; end
        dist_data = freq_IgM(:, t_idx);
        
        plot_2d_fcs(ax, values_kd, dist_data, t_val, c_code, 'IgM');
    end
    
    % --- Row 3: IgG 流式图 ---
    for i = 1:3
        % 计算位置：占据第3行的位置
        tile_idx = 12 + (i-1)*2 + 1; % 13, 15, 17
        ax = nexttile(tile_idx, [1 2]);
        
        t_val = target_times(i);
        c_code = colors{i};
        
        % 提取对应时刻数据
        t_idx = find(time_points == t_val, 1);
        if isempty(t_idx), t_idx=1; end
        dist_data = freq_IgG(:, t_idx);
        
        plot_2d_fcs(ax, values_kd, dist_data, t_val, c_code, 'IgG');
    end

    sgtitle('Spatiotemporal Dynamics of Immune Memory', 'FontSize', 16, 'FontWeight', 'bold');
end

%% --- 核心函数1：绘制带标注的 3D 曲面图 ---
function plot_3d_surface(ax, t_vec, kd_vec, freq_mat, mark_times, colors, label)
    axes(ax);
    [X, Y] = meshgrid(t_vec, kd_vec);
    
    % 绘制曲面
    s = surf(X, Y, freq_mat, 'FaceAlpha', 0.85);
    s.EdgeColor = 'none';
    s.FaceColor = 'interp';
    hold on;
    
    % 设置 Colormap
    if strcmp(label, 'IgM')
        colormap(ax, 'parula');
    else
        colormap(ax, 'jet');
    end
    
    % *** 关键修改：反转 Y 轴方向 ***
    ylim([-31, -13]);     % 设定数据范围
    set(gca, 'YDir', 'reverse'); % 让 -13 (数值较大) 显示在坐标轴的前端/上方
    
    % 绘制时间切片线 (Visual Links)
    for k = 1:length(mark_times)
        t_now = mark_times(k);
        col = colors{k};
        
        t_idx = find(t_vec == t_now, 1);
        if ~isempty(t_idx)
            z_line = freq_mat(:, t_idx);
            % 抬高一点以防遮挡
            plot3(repmat(t_now, size(kd_vec)), kd_vec, z_line + 0.01, ...
                'Color', col, 'LineWidth', 2.5);
            
            % 标注文本
            % 由于轴反转，文字位置可能需要微调
            if max(z_line) > 0.05 || t_now == 0
               text(t_now, -13.5, max(z_line)+0.15, sprintf('T=%d', t_now), ...
                    'Color', col, 'FontWeight', 'bold', 'FontSize', 9, ...
                    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
            end
        end
    end
    
    % 装饰
    xlabel('Time (days)'); 
    ylabel('Affinity (log K_d)'); 
    zlabel('Frequency');
    title([label ' Evolution'], 'FontSize', 12, 'FontWeight', 'bold');
    
    % 调整视角：这种视角配合 YDir reverse 能看清 -13 端的峰值
    view(-40, 30); 
    
    grid on; box on;
    camlight; lighting gouraud;
end

%% --- 核心函数2：绘制流式细胞散点图 ---
function plot_2d_fcs(ax, kd_vec, prob_dist, t_val, color_code, label)
    axes(ax);
    
    % 伪粒子生成
    num_particles = 3000;
    abundance = sum(prob_dist);
    % 根据丰度动态调整点数，模拟真实信号
    n_dots = round(num_particles * (abundance / 4)); 
    if n_dots > 3000, n_dots = 3000; end
    if n_dots < 10 && t_val ~= 0, n_dots = 10; end % 少许背景
    
    x_dots = []; y_dots = []; dens = [];
    
    if n_dots > 0
        cdf = cumsum(prob_dist) / sum(prob_dist);
        r = rand(n_dots, 1);
        x_base = zeros(n_dots, 1);
        
        bin_w = kd_vec(2) - kd_vec(1);
        
        for i=1:n_dots
            idx = find(cdf >= r(i), 1);
            if isempty(idx), idx = length(kd_vec); end
            % Jitter
            x_base(i) = kd_vec(idx) + (rand - 0.5) * bin_w;
        end
        x_dots = x_base;
        y_dots = 400 + randn(n_dots, 1) * 150; % SSC 模拟
        y_dots(y_dots < 0) = 0; y_dots(y_dots > 1000) = 1000;
        dens = x_dots; 
    end
    
    scatter(x_dots, y_dots, 8, dens, 'filled', 'MarkerFaceAlpha', 0.6);
    
    if strcmp(label, 'IgM')
        colormap(ax, 'parula');
    else
        colormap(ax, 'jet');
    end
    % 保持流式图的X轴按照数值习惯 (左小右大)
    % 这样可以看出三维图虽然翻转了展示视角，但二维分析依然符合数轴逻辑
    xlim([-31.5, -12.5]); 
    ylim([0, 1000]);
    caxis([-31 -13]);
    
    % 视觉美化
    title(sprintf('%s at T=%d', label, t_val), 'Color', color_code, 'FontWeight', 'bold');
    
    if t_val == 0 
        ylabel('SSC');
    else
        set(gca, 'YTickLabel', []);
    end
    xlabel('log K_d');
    
    if isempty(x_dots)
        text(-22, 500, 'Not Detected', 'Color', [0.6 0.6 0.6], 'HorizontalAlignment', 'center');
    end
    
    box on; grid on;
    set(gca, 'LineWidth', 1.1, 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2]);
end

%% --- 数据生成与处理部分 (保持不变) ---

function [t, kd, d_m, d_g] = run_simulation()
    % 参数设置
    M_1 = 1e15; M_2 = 4e18;
    G_1 = 1e15; G_2 = 4e19; E = 1e18;
    C_1 = 1e13; C_2 = 4e16; C_3 = 1e13; C_4 = 4e17;
    P_M = 5e7;  P_G = 1e8;

    % 初始分布 prob_A
    mu = 9; sigma = 0.8;
    prob_A = zeros(1,10);
    prob_A(1) = normcdf(5, mu, sigma);
    for i = 2:5
        prob_A(i) = normcdf(5+i-1, mu, sigma) - normcdf(5+i-2, mu, sigma);
    end
    for i = 6:10, prob_A(i) = prob_A(11-i); end
    
    AA = zeros(100,1);
    for i=1:10, for j=1:10, AA(10*(i-1)+j) = prob_A(i)*prob_A(j); end, end

    % x0 Initialization
    x0 = zeros(1402,1);
    for i=1:100
        term = AA(i);
        x0(i) = M_1 * term; x0(i+100) = M_2 * term;
        x0(i+200) = G_1 * term; x0(i+300) = G_2 * term;
        x0(i+400) = C_1 * term; x0(i+500) = C_2 * term;
        x0(i+600) = C_3 * term; x0(i+700) = C_4 * term;
        x0(i+1200)= P_M * term; x0(i+1300)= P_G * term;
    end
    x0(801:1200) = 0; x0(1401) = E; x0(1402) = 10;
    
para(1) = 1e-20; % environmental antigen kon
para(2) = 0.5;% environmental antigen koff
para(3) = 2.2e17+1e13; % replenish constant pi 1 
para(4) = 0.5e13;% replenish constant pi 2
para(5) = 0.01; % decay constant of BCR IgM
para(6) = 0.005;% decay constant of BCR IgG
para(7) = 1; % k2  feedback constant of enviromental antigen-antibody complex
para(8) = 5e-7;% k2' feedback constant on PC cell regeneration
para(9) = 0.1; % decay constant of plasma Cell IgM
para(10) = 0.05;%% decay constant of plasma Cell IgG
para(11) = 4.4e9;% production constant of IgM
para(12) = 1.2e10;% production constant of IgG

para(13) = 0.05;% decay constant of IgM
para(14) = 0.025;% decay constant of IgG
para(15) = 2;% amplification constant of virus antigen
para(16) = 0.02;% transformation constant from IgM to IgG memory cell
para(17) = 0.1;% maximal production percentage of plasma cell 
para(18) = 0.5;% decay constant of complex
para(19) = 1.2; % virus replication constant
para(20) = 1e5;



para_new(1) = 1e-22; 
para_new(2) = 1e-21;
para_new(3) = 1e-20; 
para_new(4) = 1e-19;
para_new(5) = 1e-18;
para_new(6) = 1e-17;
para_new(7) = 1e-16; 
para_new(8) = 1e-15;
para_new(9) = 1e-14; 
para_new(10) = 1e-13;

para_new_1(1) = 1e0; 
para_new_1(2) = 1e1;
para_new_1(3) = 1e2; 
para_new_1(4) = 1e3;
para_new_1(5) = 1e4; 
para_new_1(6) = 1e5;
para_new_1(7) = 1e6; 
para_new_1(8) = 1e7;
para_new_1(9) = 1e8; 
para_new_1(10) = 1e9;
    
    [t, y] = ode15s(@pathway_model_many_antibody_immune_include_plasma, [0 1000], x0, [], para, para_new, para_new_1, AA);
    
    t_interp = 0:5:1000;
    raw_m = interp1(t, y(:,101:200), t_interp);
    raw_g = interp1(t, y(:,301:400), t_interp);
    
    % Map to KD bins
    d_m = zeros(19, length(t_interp));
    d_g = zeros(19, length(t_interp));
    
    for i=1:100
        % 原始逻辑: fix((i-1)/10) - mod(i-1,10) + 10
        row_idx = fix((i - 1)/10) - mod(i - 1, 10) + 10;
        if row_idx >=1 && row_idx <= 19
            d_m(row_idx, :) = d_m(row_idx, :) + raw_m(:, i)';
            d_g(row_idx, :) = d_g(row_idx, :) + raw_g(:, i)';
        end
    end
    t = t_interp;
    kd = -31:1:-13;
    
    d_m = d_m / max(d_m(:));
    d_g = d_g / max(d_g(:));
end

function [t, kd, d_m, d_g] = generate_mock_data()
    t = 0:5:1000;
    kd = -31:1:-13; 
    d_m = zeros(19, length(t));
    d_g = zeros(19, length(t));
    
    % Mock IgM: 早期峰值偏向 -25，但也覆盖 -13
    for i=1:19
        k_val = kd(i);
        af_factor = exp(-(k_val - (-22))^2 / 15);
        tm_factor = exp(-(t-50).^2 / 3000);
        d_m(i,:) = af_factor * tm_factor;
    end
    
    % Mock IgG: 
    for i=1:19
        k_val = kd(i);
        af_factor = exp(-(k_val - (-18))^2 / 10); 
        tm_factor = 1 ./ (1 + exp(-0.03*(t-120)));
        d_g(i,:) = af_factor * tm_factor;
    end
end