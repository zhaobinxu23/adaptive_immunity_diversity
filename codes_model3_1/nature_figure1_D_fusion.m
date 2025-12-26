function nature_figure1_D_fusion()
    % 清除环境
    clc; clear; close all;

    % ==========================================
    % 1. 参数与初始条件设置 (完全保持原样)
    % ==========================================
    para = zeros(11,1);
    para(1) = 1e-7;   % k1 drug binding
    para(2) = 1e-14;  % k-1 drug dissociation
    para(3) = 1e-2;   % k2 drug decay
    para(4) = 0.5;    % k3 complex decay
    para(5) = 1;      % k4 virus replication
    para(6) = 1e-7;   % k5 antibody binding
    para(7) = 1e-14;  % k-5 antibody dissociation
    para(8) = 10;     % k6 antibody replenish
    para(9) = 1e-2;   % k7 antibody decay
    para(10) = 2;     % k8 feedback
    para(11) = 1e-1;  % p

    % 初始状态 Phase 1
    x0 = [0; 1; 0; 1000; 0; 0]; 

    % ==========================================
    % 2. 模拟运行
    % ==========================================
    
    % --- Phase 1: 0-20 天 (无药) ---
    options = odeset('RelTol',1e-6, 'AbsTol',1e-9);
    [t1, y1] = ode15s(@pathway_model_with_drug_inhib, [0 20], x0, options, para);

    % --- Phase 2: 20-200 天 (给药) ---
    % 继承上一阶段的状态
    x_new = y1(end, :)';
    
    % !!! 关键干预：注入药物 !!!
    x_new(1) = 7e8; 
    
    [t2, y2] = ode15s(@pathway_model_with_drug_inhib, [20 200], x_new, options, para);

    % --- 合并数据 ---
    t_total = [t1; t2];
    y_drug  = [y1(:,1); y2(:,1)];
    y_virus = [y1(:,2); y2(:,2)];
    y_ab    = [y1(:,4); y2(:,4)];

    % ==========================================
    % 3. Nature 风格绘图
    % ==========================================
    
    % 定义配色 (RGB格式)
    c_vir  = [230, 75, 53] / 255;    % 红色 (Virus)
    c_ab   = [0, 160, 135] / 255;    % 蓝绿色 (Antibody)
    c_drug = [60, 84, 136] / 255;    % 深蓝色 (Drug)
    
    % 创建高分辨率画布
    figure('Color', 'w', 'Units', 'pixels', 'Position', [200, 200, 900, 600]);
    
    % --- 激活双坐标轴 ---
    yyaxis right  % 先画右轴 (背景层：药物)
    ax_right = gca;
    
    % 使用 area 绘制半透明的药物浓度背景，非常有高级感
    area(t_total, y_drug, 'FaceColor', c_drug, 'EdgeColor', 'none', 'FaceAlpha', 0.15);
    hold on;
    % 叠加一条药物轮廓虚线
    p_drug = plot(t_total, y_drug, '--', 'Color', c_drug, 'LineWidth', 1.5);
    
    % 右轴设置
    ylabel('Drug Concentration (nM)', 'FontSize', 12, 'FontWeight', 'bold');
    set(ax_right, 'YColor', c_drug, 'YLim', [0, 8e8]); % 线性坐标展示药物消减
    
    
    % --- 切换到左轴 (核心数据层) ---
    yyaxis left
    ax_left = gca;
    
    % 绘制病毒 (红色加粗实线)
    p_vir = plot(t_total, y_virus, '-', 'Color', c_vir, 'LineWidth', 3);
    hold on;
    % 绘制抗体 (绿色点线)
    p_ab = plot(t_total, y_ab, ':', 'Color', c_ab, 'LineWidth', 2.5);
    
    % 左轴设置 (关键：使用对数坐标 Log Scale)
    set(ax_left, 'YScale', 'log', 'YColor', [0.2 0.2 0.2]); % 轴色设为深灰
    ylabel('Viral Load / Antibody Titers (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold');
    ylim([1e-1, 1e8]); % 设置合理的显示范围
    xlim([0, 200]);

    % ==========================================
    % 4. 标注与美化 (Annotations)
    % ==========================================
    
    % 添加一条淡灰色的竖线标记给药时间
    xline(20, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    text(22, 2e7, 'Drug Injection (Day 20)', 'Color', [0.4 0.4 0.4], 'FontSize', 10, 'FontWeight', 'bold');

    % 标注病毒反弹 (Rebound)
    % 自动寻找第二个峰值
    [max_val_ph2, idx_ph2] = max(y2(:,2));
    t_peak = t2(idx_ph2);
    
    if max_val_ph2 > 1e2 % 只有存在明显反弹才标注
        % 画一个醒目的箭头指向峰值
        text(t_peak, max_val_ph2*2, {'Viral Rebound'}, ...
            'Color', c_vir, 'HorizontalAlignment', 'center', ...
            'FontSize', 11, 'FontWeight', 'bold');
        plot(t_peak, max_val_ph2, 'o', 'Color', c_vir, 'MarkerSize', 8, 'LineWidth', 2);
    end

    % 标注抗体水平低的原因 (Sub-threshold)
    % 在抗体低谷处添加简短说明
    text(80, 2e0, {'\it suppressed antigen', '\it limits antibody response'}, ...
         'Color', c_ab, 'FontSize', 9, 'HorizontalAlignment', 'center');

    % 通用修饰
    xlabel('Time (Days)', 'FontSize', 12, 'FontWeight', 'bold');
    title('Dynamics of Viral Rebound following Inhibitor Therapy', 'FontSize', 14, 'FontName', 'Arial');
    
    % 图例优化
    legend([p_vir, p_ab, p_drug], {'Viral Load', 'Endogenous Antibody', 'Drug Concentration'}, ...
           'Location', 'northeast', 'Box', 'off', 'FontSize', 10);
       
    % 去掉丑陋的边框，现代风格
    set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11, 'FontName', 'Arial');
    grid off;

    hold off;
end

% ==========================================
% 5. 你的原始方程函数 (稍微整理格式)
% ==========================================
function F = pathway_model_with_drug_inhib(t,y,para)
    % 确保非负
    y = max(0,y);
    
    F = zeros(6,1); 
    
    % 变量赋值
    D = y(1); V = y(2); DV = y(3); 
    Ab = y(4); AbV = y(5); DVAb = y(6); % 修正变量名对应关系，保持你原始逻辑
    
    % 参数赋值
    k1=para(1); k_1=para(2); k2=para(3); k3=para(4); k4_repl=para(5);
    k5_bind=para(6); k_5_off=para(7); k6_base=para(8); k7_decay=para(9); 
    k8_feed=para(10); p_val=para(11);

    % 方程组 (完全对应你给出的 F(1,1) 到 F(6,1))
    
    % Drug
    F(1) = -k1*D*V + k_1*DV - k2*D;
    
    % Virus
    F(2) = -k1*D*V + k_1*DV + k4_repl*p_val*DV + k4_repl*V - k5_bind*V*Ab + k_5_off*AbV;
    
    % DV Complex
    F(3) = k1*D*V - k_1*DV - k5_bind*DV*Ab + k_5_off*DVAb;
    
    % Antibody
    F(4) = -k5_bind*(V+DV)*Ab + k_5_off*(AbV+DVAb) + k6_base - k7_decay*Ab + k8_feed*(AbV+DVAb);
    
    % AbV Complex
    F(5) = k5_bind*V*Ab - k_5_off*AbV - k3*AbV;
    
    % DVAb Complex
    F(6) = k5_bind*DV*Ab - k_5_off*DVAb - k3*DVAb;
end