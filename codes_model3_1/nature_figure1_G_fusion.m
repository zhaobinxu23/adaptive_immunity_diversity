function nature_figure1_G_fusion
    % -------------------------------------------------------------------------
    % Combined Immunological Dynamics Simulation
    % Developed for Nature-standard visualization
    % -------------------------------------------------------------------------
    
    clc; clear; close all;

    % --- Global Plotting Settings ---
    fig = figure('Units', 'pixels', 'Position', [100, 100, 1000, 800]);
    set(fig, 'Color', 'w'); % White background
    
    % Define Colors (Nature-style compatible)
    color_strong_virus = [0.8500, 0.3250, 0.0980]; % Orange/Red
    color_weak_virus   = [0.0000, 0.4470, 0.7410]; % Dark Blue
    
    lw = 2.5;      % LineWidth
    fs = 12;       % FontSize
    font_name = 'Arial';
    
    % Setup Tiled Layout (2 Rows, 2 Columns for key comparisons)
    t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    
    %% --- SCENARIO 1: Strong Immunity (k2 = 3) ---
    % Simulate Strong Virus
    [t1_s, y1_s] = run_model(3, 0.55); 
    % Simulate Weak Virus
    [t1_w, y1_w] = run_model(3, 0.50);
    
    % Panel A: Antibody Dynamics (Strong Immunity)
    nexttile;
    p1 = plot(t1_s, y1_s(:,1), 'Color', color_strong_virus, 'LineWidth', lw); hold on;
    p2 = plot(t1_w, y1_w(:,1), 'Color', color_weak_virus, 'LineWidth', lw, 'LineStyle', '--');
    title('A. Antibody (Strong Immunity)', 'FontWeight', 'bold');
    ylabel('Concentration', 'FontName', font_name);
    grid on; box on;
    set(gca, 'FontSize', fs, 'FontName', font_name, 'LineWidth', 1.2);
    legend([p1, p2], {'Strong Virus (k3=0.55)', 'Weak Virus (k3=0.5)'}, ...
           'Location', 'Best', 'Box', 'off');

    % Panel B: Virus Dynamics (Strong Immunity)
    nexttile;
    plot(t1_s, y1_s(:,5), 'Color', color_strong_virus, 'LineWidth', lw); hold on;
    plot(t1_w, y1_w(:,5), 'Color', color_weak_virus, 'LineWidth', lw, 'LineStyle', '--');
    title('B. Virus Load (Strong Immunity)', 'FontWeight', 'bold');
    grid on; box on;
    set(gca, 'FontSize', fs, 'FontName', font_name, 'LineWidth', 1.2);

    
    %% --- SCENARIO 2: Weak Immunity (k2 = 2) ---
    % Simulate Strong Virus
    [t2_s, y2_s] = run_model(2, 0.55); 
    % Simulate Weak Virus
    [t2_w, y2_w] = run_model(2, 0.50); 

    % Panel C: Antibody Dynamics (Weak Immunity)
    nexttile;
    plot(t2_s, y2_s(:,1), 'Color', color_strong_virus, 'LineWidth', lw); hold on;
    plot(t2_w, y2_w(:,1), 'Color', color_weak_virus, 'LineWidth', lw, 'LineStyle', '--');
    title('C. Antibody (Weak Immunity)', 'FontWeight', 'bold');
    xlabel('Time (hours/days)', 'FontName', font_name);
    ylabel('Concentration', 'FontName', font_name);
    grid on; box on;
    set(gca, 'FontSize', fs, 'FontName', font_name, 'LineWidth', 1.2);
    
    % Panel D: Virus Dynamics (Weak Immunity)
    nexttile;
    plot(t2_s, y2_s(:,5), 'Color', color_strong_virus, 'LineWidth', lw); hold on;
    plot(t2_w, y2_w(:,5), 'Color', color_weak_virus, 'LineWidth', lw, 'LineStyle', '--');
    title('D. Virus Load (Weak Immunity)', 'FontWeight', 'bold');
    xlabel('Time (hours/days)', 'FontName', font_name);
    grid on; box on;
    set(gca, 'FontSize', fs, 'FontName', font_name, 'LineWidth', 1.2);
    
    % Add a main title for the whole figure
    title(t, 'Comparison of Immune Dynamics under Variable Viral Loads', ...
        'FontSize', 16, 'FontName', font_name, 'FontWeight', 'bold');

end 

%% --- Underlying Helper Functions ---

function [t, y] = run_model(k2_val, k3_val)
    % Setup Initial Conditions
    x0 = zeros(6,1);
    x0(1) = 0;   % A (Antibody)
    x0(2) = 500; % B (Naive B cells)
    x0(3) = 0;   % C1
    x0(4) = 0;   % C2
    x0(5) = 10;  % V (Virus initial load)
    x0(6) = 0;   % ASC (Plasma cells)

    % Base Parameters
    para(1) = 1e-7;   % k1 (Binding)
    para(2) = 1e-14;  % k-1 (Dissociation)
    para(3) = k2_val; % k2 (Immune Response Strength) <--- Variable
    para(4) = 1e-3;   % alpha
    para(5) = 0.005;  % d1 (Decay B)
    para(6) = 0.02;   % d2 (Decay A)
    para(7) = 0.5;    % d3 (Decay C1)
    para(8) = 0.5;    % d4 (Decay C2)
    para(9) = k3_val; % k3 (Virus Replication) <--- Variable
    para(10) = 1e-7;  % k4
    para(11) = 0.01;  % decay ratio of ASC cell
    para(12) = 1e5;   % Antibody secretion rate
    para(13) = 0.5;   % Masking effect / Feedback

    % Run ODE
    [t, y] = ode15s(@(t,y) pathway_model_eq(t,y,para), [0 150], x0);
end

function F = pathway_model_eq(t, y, para)
    F = zeros(6,1);
    
    % Map variables for clarity
    A = y(1); B = y(2); C1 = y(3); C2 = y(4); V = y(5); ASC = y(6);
    
    % Map Parameters
    k1 = para(1); k_m1 = para(2); k2 = para(3); alpha = para(4);
    d1 = para(5); d2 = para(6); d3 = para(7); d4 = para(8);
    k3 = para(9); k4 = para(10); d_ASC = para(11); k_sec = para(12); 
    mask = para(13);

    % Equations 
    % 1. Antibody (A) 
    % Note: Reconstructed standard eq based on logic (Production - Binding + Dissoc - Decay)
    F(1) = k_sec * ASC - k1 * A * V + k_m1 * C1 - d2 * A; 
    
    % 2. B Cells (B) (From Context)
    % Term: masking effect inhibits expansion
    F(2) = k2 * (1 - mask * C2 / (C2 + B + 1e-10)) * C2 - k4 * B * V + k_m1 * C2 - d1 * B;
    
    % 3. Complex 1 (C1: Antibody-Virus)
    F(3) = k1 * A * V - k_m1 * C1 - d3 * C1;
    
    % 4. Complex 2 (C2: BCR-Virus)
    F(4) = k4 * B * V - k_m1 * C2 - d4 * C2;
    
    % 5. Virus (V)
    F(5) = k3 * V - k1 * A * V - k4 * B * V + k_m1 * (C1 + C2);
    
    % 6. ASC Cells
    % Context Equation: para(4)*para(3)*para(13)*y(4)/(y(4)+y(2))*y(4)-para(11)*y(6)
    F(6) = alpha * k2 * mask * (C2 / (C2 + B + 1e-10)) * C2 - d_ASC * ASC;
end