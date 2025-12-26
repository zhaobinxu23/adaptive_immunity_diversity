function nature_figure1_C_fusion()
    % 清除环境
    clear; clc; close all;

    % ==========================================
    % 1. 参数设置与模型模拟
    % ==========================================
    % 对应文档中的 para(1) - para(8)
    % p1: Endo-Ab binding, p6: mAb binding, etc.
    para = [1e-5, 1e-14, 1.5, 0.5, 0.5, 1e-5, 1e-14, 0.05];

    % --- Phase 1: 0 - 28 天 (自然感染期) ---
    x0 = [10; 0; 1; 0; 0]; % 初始条件 [EndoAb; mAb; Virus; Cplx1; Cplx2]
    tspan1 = linspace(0, 28, 200);
    
    % 使用 ode15s 求解刚性方程
    [t1, y1] = ode15s(@(t,y) pathway_model(t,y,para), tspan1, x0);

    % --- Phase 2: 28 - 150 天 (治疗介入后) ---
    % 提取第28天的状态作为新的初始值
    x_intervention = y1(end, :)';
    % 注入单克隆抗体 (mAb = 1.5e6)
    x_intervention(2) = 1.5e6;

    tspan2 = linspace(28, 150, 1000);
    [t2, y2] = ode15s(@(t,y) pathway_model(t,y,para), tspan2, x_intervention);

    % 合并数据
    t_total = [t1; t2];
    y_endo_ab = [y1(:, 1); y2(:, 1)];
    y_mab = [y1(:, 2); y2(:, 2)];
    y_virus = [y1(:, 3); y2(:, 3)];

    % ==========================================
    % 2. Nature 风格配色与绘图设置
    % ==========================================
    % 定义 NPG (Nature Publishing Group) 风格颜色 (RGB 归一化到 0-1)
    color_virus = [230, 75, 53] / 255;   % 红色
    color_mab   = [77, 187, 213] / 255;  % 蓝色
    color_endo  = [0, 160, 135] / 255;   % 绿色
    color_grey  = [0.5, 0.5, 0.5];

    % 创建图形窗口
    figure('Color', 'w', 'Position', [100, 100, 900, 600]);
    ax = gca;
    hold on;

    % 设置对数坐标轴范围 (避免 log(0) 和适应 quantitive mismatch)
    y_min_limit = 1e-2; 
    y_max_limit = 1e7;
    ylim([y_min_limit, y_max_limit]);
    xlim([0, 150]);

    % --- A. 绘制背景机制区域 (Kinetic Decoupling Gap) ---
    % 在 mAb 衰减且内源抗体未升起的时间段画灰色背景 (约 60-100天)
    % 注意：patch 的 Y 坐标要覆盖整个 Y 轴范围
    fill_x = [60, 100, 100, 60];
    fill_y = [y_min_limit, y_min_limit, y_max_limit, y_max_limit];
    patch(fill_x, fill_y, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
    
    % 添加机制说明文字
    text(80, 2e-2, {'Mechanism: Kinetic Decoupling', '\it(mAb decays before Endo-Ab maturation)'}, ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'Color', [0.4 0.4 0.4], 'Interpreter', 'tex');

    % --- B. 绘制主要曲线 (使用半对数坐标逻辑) ---
    % MATLAB 的 semilogy 实际上是 plot 但改了坐标轴属性，这里我们手动控制
    
    % 1. 病毒载量 (红实线)
    p1 = plot(t_total, y_virus, 'Color', color_virus, 'LineWidth', 2.5);
    
    % 2. 单克隆抗体 (蓝虚线)
    p2 = plot(t_total, y_mab, '--', 'Color', color_mab, 'LineWidth', 2);
    
    % 3. 内源性抗体 (绿点线)
    p3 = plot(t_total, y_endo_ab, ':', 'Color', color_endo, 'LineWidth', 2);

    % 设置为对数坐标
    set(ax, 'YScale', 'log');

    % --- C. 关键事件标注 (Storytelling) ---

    % 1. 注射时刻 (Day 28)
    xline(28, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 1.5, 'Alpha', 0.7);
    text(29, 2e6, {'\bf mAb Injection', '(Day 28)'}, 'VerticalAlignment', 'top', 'Color', [0.2 0.2 0.2], 'FontSize', 10);

    % 2. 病毒反弹 (Viral Rebound) 标注
    % 找到 Phase 2 中的峰值
    [peak_val, peak_idx] = max(y2(:, 3));
    peak_time = t2(peak_idx);
    
    % 只有当反弹明显时才标注
    if peak_time > 60 && peak_val > 1
        % 绘制箭头标注
        annotation('textarrow', [0.65, 0.58], [0.6, 0.5], ... % 归一化坐标，需根据图手动微调或使用数据坐标转换
            'String', {'Viral', 'Rebound'}, 'Color', color_virus, 'FontSize', 11, 'FontWeight', 'bold', 'TextColor', color_virus);
        
        % 在峰值处画一个虚线圆圈强调
        plot(peak_time, peak_val, 'o', 'MarkerSize', 25, 'Color', color_virus, 'LineWidth', 1, 'LineStyle', '--');
    end

    % 3. 阈值线 (只是示意)
    yline(10, '-.', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    text(145, 12, 'Neutralization Threshold', 'HorizontalAlignment', 'right', 'FontSize', 9, 'Color', [0.5 0.5 0.5]);

    % ==========================================
    % 3. 图表排版与美化
    % ==========================================
    
    % 坐标轴标签
    xlabel('Time (Days)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial');
    ylabel('Concentration / Viral Load (Log Scale)', 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial');
    
    % 图例设置
    legend([p1, p2, p3], {'Viral Load', 'Monoclonal Ab (mAb)', 'Endogenous Ab'}, ...
        'Location', 'northeast', 'Box', 'off', 'FontSize', 10, 'FontName', 'Arial');
    
    % 调整坐标轴外观 (Nature 风格：去掉上边和右边的框线)
    set(ax, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1.2, 'FontSize', 11, 'FontName', 'Arial');
    
    % 标题
    title('Kinetic Dynamics of Viral Rebound under mAb Therapy', 'FontSize', 14, 'FontName', 'Arial', 'FontWeight', 'bold');

    hold off;
end

% ==========================================
% 4. ODE 模型定义函数
% ==========================================
function dy = pathway_model(t, y, para)
    % y(1): Endo Ab
    % y(2): Injected mAb
    % y(3): Virus
    % y(4): Endo Ab-Virus Complex
    % y(5): mAb-Virus Complex

    % 确保非负
    y = max(0, y);

    % 参数解包
    p1 = para(1); p2 = para(2); p3 = para(3); p4 = para(4);
    p5 = para(5); p6 = para(6); p7 = para(7); p8 = para(8);

    dy = zeros(5,1);
    
    % d(Endo_Ab)/dt
    dy(1) = -p1*y(1)*y(3) + p2*y(4) + p3*y(4);
    
    % d(mAb)/dt
    dy(2) = -p6*y(2)*y(3) + p7*y(5) - p8*y(2);
    
    % d(Virus)/dt
    dy(3) = -p1*y(1)*y(3) + p2*y(4) - p6*y(2)*y(3) + p7*y(5) + p5*y(3);
    
    % d(Endo_Cplx)/dt
    dy(4) = p1*y(1)*y(3) - p2*y(4) - p4*y(4);
    
    % d(mAb_Cplx)/dt
    dy(5) = p6*y(2)*y(3) - p7*y(5) - p4*y(5);
end