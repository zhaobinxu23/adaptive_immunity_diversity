function nature_figure1_F_fusion()
  
   
    % -------------------------------------------------------------------------
    % Figure 1: Effect of Vaccination Interval and Dosage on Antibody Response
    % Optimized for Nature-style visualization
    % -------------------------------------------------------------------------
    
    clc; clear; close all;

    % --- Figure Setup ---
    fig = figure('Units', 'pixels', 'Position', [100, 100, 1200, 400]);
    set(fig, 'Color', 'w'); 
    
    % Create Tiled Layout (1 row, 3 columns)
    t_layout = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(t_layout, 'Impact of Vaccination Interval and Antigen Dosage on Antibody Response', ...
        'FontSize', 16, 'FontName', 'Arial', 'FontWeight', 'bold');

    % Common Parameters
    % k1, k_m1, k2, alpha, d1, d2, d3, d4, k3, k4, d_ASC, k_sec, mask
    para = [1e-7, 1e-14, 2, 1e-3, 0.005, 0.02, 0.5, 0.5, -0.02, 1e-7, 0.01, 1e5, 0.5];
    
    % Initial State: [A, B, C1, C2, V, ASC]
    x0_init = [0; 500; 0; 0; 1e8; 0]; 

    %% --- Scenario 1: Long Interval (500 time units) ---
    % Step 1: First Shot
    [t1a, y1a] = ode15s(@(t,y) model_eq(t,y,para), [0 500], x0_init);
    
    % Step 2: Second Shot (Standard Dose 1e8)
    x0_2nd = y1a(end, :)';
    x0_2nd(5) = 1e8; % Re-inject standard dose
    [t1b, y1b] = ode15s(@(t,y) model_eq(t,y,para), [500 1000], x0_2nd);
    
    % Plot Panel A
    nexttile;
    plot_scenario([t1a; t1b], [y1a; y1b], 'A. Long Interval (500 t)', 'Standard Dose (1e8)');
    
    
    %% --- Scenario 2: Short Interval (100 time units) ---
    % Step 1: First Shot
    [t2a, y2a] = ode15s(@(t,y) model_eq(t,y,para), [0 100], x0_init);
    
    % Step 2: Second Shot (Standard Dose 1e8) - The "Interference" Scenario
    x0_2nd = y2a(end, :)';
    x0_2nd(5) = 1e8; % Re-inject standard dose
    [t2b, y2b] = ode15s(@(t,y) model_eq(t,y,para), [100 600], x0_2nd);
    
    % Plot Panel B
    nexttile;
    plot_scenario([t2a; t2b], [y2a; y2b], 'B. Short Interval (100 t)', 'Standard Dose (1e8)');
    % Highlight specific area
    xline(100, '--k', 'Alpha 0.3'); 
    text(120, max(y2b(:,1))*0.5, {'Interference:', 'Weak Boost'}, 'Color', 'red', 'FontSize', 10);

    %% --- Scenario 3: Short Interval + High Dose (The Solution) ---
    % Step 1: First Shot
    [t3a, y3a] = ode15s(@(t,y) model_eq(t,y,para), [0 100], x0_init);
    
    % Step 2: Second Shot (High Dose 1e10) - The "Breakthrough" Scenario
    x0_2nd = y3a(end, :)';
    x0_2nd(5) = 1e10; % Re-inject HIGH dose
    [t3b, y3b] = ode15s(@(t,y) model_eq(t,y,para), [100 600], x0_2nd);
    
    % Plot Panel C
    nexttile;
    plot_scenario([t3a; t3b], [y3a; y3b], 'C. Short Interval + High Dose', 'High Dose (1e10)');
    
end

%% --- Helper Functions ---

function plot_scenario(t, y, list_title, list_legend)
    % Consistent plotting style
    plot(t, y(:,1), 'LineWidth', 2.5, 'Color', [0.2, 0.4, 0.7]); % Blue line for Antibody
    grid on; box on;
    title(list_title, 'FontSize', 12, 'FontWeight', 'bold');
    xlabel('Time', 'FontSize', 11);
    ylabel('Antibody Titer', 'FontSize', 11);
    legend({list_legend}, 'Location', 'northeast', 'Box', 'off');
    xlim([0, max(t)]);
    ylim([0, inf]);
    set(gca, 'LineWidth', 1.2, 'FontSize', 11, 'FontName', 'Arial');
end

function F = model_eq(~, y, para)
    % Unpack Parameters
    k1 = para(1); k_m1 = para(2); k2 = para(3); alpha = para(4);
    d1 = para(5); d2 = para(6); d3 = para(7); d4 = para(8);
    k3 = para(9); k4 = para(10); d_ASC = para(11); k_sec = para(12);
    mask = para(13);

    % Variables: A(1), B(2), C1(3), C2(4), V(5), ASC(6)
    F = zeros(6,1);
    
    % Non-negative constraint logic (optional within solver, but implemented as max(0,y) in calling not ideal inside ODE)
    % Assuming standard dynamics here.
    
    % 1. Antibody (A)
    F(1) = k_sec * y(6) - k1 * y(1) * y(5) + k_m1 * y(3) - d2 * y(1);
    
    % 2. Naive B cells (B)
    % Inhibition term: mask * C2 / (C2 + B)
    inhibition = 1 - mask * y(4) / (y(4) + y(2) + 1e-20);
    F(2) = k2 * inhibition * y(4) - k4 * y(2) * y(5) + k_m1 * y(4) - d1 * y(2);
    
    % 3. Complex 1 (Ab-V)
    F(3) = k1 * y(1) * y(5) - k_m1 * y(3) - d3 * y(3);
    
    % 4. Complex 2 (BCR-V)
    F(4) = k4 * y(2) * y(5) - k_m1 * y(4) - d4 * y(4);
    
    % 5. Virus/Antigen (V)
    F(5) = k3 * y(5) - k1 * y(1) * y(5) - k4 * y(2) * y(5) + k_m1 * (y(3) + y(4));
    
    % 6. ASC
    % Activation term
    F(6) = alpha * k2 * mask * (y(4) / (y(4) + y(2) + 1e-20)) * y(4) - d_ASC * y(6);
end