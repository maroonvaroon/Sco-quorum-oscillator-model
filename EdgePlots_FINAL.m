close all; clear;
tic
%  ===================
%  SPEED / RESOLUTION
%  ===================
dx_base = 1; dt_base = 0.1;
dx_ref  = 0.125;
T_ref = 49000;

tau_cal = 0.7057;

dx = 0.25;

s  = dx / dx_base;
dt = s^2 * dt_base;

tFinal_ref = T_ref * ( (dx_ref/dx_base)^2 * dt_base );
T = max(1, round(tFinal_ref / dt));

PLOT_EVERY = 20;
KYMO_EVERY = 5;

USE_LAPLACIAN_9 = true;

fprintf('dx=%.4f, dt=%.6f, T=%d (keeps ~same physical horizon)\n', dx, dt, T);

%  =======================
%  DOMAIN DIMENSIONS/SIZE
%  =======================
Lx = 60; Ly = 60;
Nx = max(4, round(Lx / dx)); Ny = max(4, round(Ly / dx));

%  ==========================
%  MODEL PARAMETERS
%  ==========================

% Sharpness of repression

n = 10;

% QS threshold

tau = 0.1;

% Copy number for bacterial growth / Growth rate

mu = 0.3;

% Sharpness of QS filter

k = 100;

% Growth Process

K_N = 2000;  % Carrying Capacity
D_0 = 10000;    % Colony Growth "Diffusion"/Spread
eps_border = 0.01; % Sharpness of the border (smaller=larger)

% Diffusion Coefficients

D_A  = 0;
D_G  = 1;
D_R  = 0;
D_R2 = 0;
D_C  = 0;

% Production Coefficients

alpha_A  = 1;
alpha_G  = 1;
alpha_R  = 1.1;
alpha_R2 = 1.8;
alpha_C  = 3;

% Degradation Coefficients

delta_A  = 1;
delta_G  = 1;
delta_R  = 0.1;
delta_R2 = 0.5;
delta_C  = 0.25;

% Repression Strengths

K_R2A = 0.4;
K_RR  = 0.5;
K_RC  = 1;
K_RR2 = 0.5;
K_R2C = 1;

%  ==========================
%  INITIAL CONDITIONS
%  ==========================
[X, Y] = meshgrid(1:Nx, 1:Ny);
x0 = Nx / 2; y0 = Ny / 2;
cx = round(x0); cy = round(y0); 

sigma_phys = 3; sigma = sigma_phys / dx;
N  = K_N * 0.01 * exp(-((X - x0).^2 + (Y - y0).^2) / (2 * sigma^2));
A = zeros(Ny, Nx); G = A; R = A; R2 = A; C = A;
A(cy, cx) = 0.1;

%  =================================
%  DYNAMIC EDGE TRACKING PARAMETERS
%  =================================
R_cells = hypot(X - cx, Y - cy) * dx;
N_edge_threshold = 0.05 * K_N; 

% Adjust ring_width below to average over a wider or narrower ring
ring_width = 0.5; 

if USE_LAPLACIAN_9
    lap = @(U) (0.25*(circshift(U,[1,1])+circshift(U,[-1,-1])+circshift(U,[1,-1])+circshift(U,[-1,1])) + ...
        0.5*(circshift(U,[1,0])+circshift(U,[-1,0])+circshift(U,[0,1])+circshift(U,[0,-1])) - 3*U) / dx^2;
else
    lap = @(U) (circshift(U,[1,0])+circshift(U,[-1,0])+circshift(U,[0,1])+circshift(U,[0,-1]) - 4*U) / dx^2;
end

%  ==========================================
%  DATA PRE-ALLOCATION
%  ==========================================
time_vector = (0:T-1) * dt*tau_cal;
N_edge  = zeros(T, 1);
A_edge  = zeros(T, 1);
G_edge  = zeros(T, 1);
R_edge  = zeros(T, 1);
R2_edge = zeros(T, 1);
C_edge  = zeros(T, 1);

%% ================
%  TIME LOOP 
%  ================
for ti = 1:T
    H = 1 ./ (1 + exp(-(N - 2 * K_N) / eps));
    D = D_0 * H;
    dNdx = (circshift(N, [0, -1]) - circshift(N, [0, 1])) / (2*dx);
    dNdy = (circshift(N, [-1, 0]) - circshift(N, [1, 0])) / (2*dx);
    dDdx = (circshift(D, [0, -1]) - circshift(D, [0, 1])) / (2*dx);
    dDdy = (circshift(D, [-1, 0]) - circshift(D, [1, 0])) / (2*dx);
    div_grad = D .* lap(N) + dDdx .* dNdx + dDdy .* dNdy;
    growth   = mu * N .* (1 - N / K_N);
    N        = N + dt * (div_grad + growth);
    N        = max(N, 0);
    N        = min(N, 1.1 * K_N);
    
    lower_threshold = tau/10;
    upper_threshold = 7.5*tau;
    H_G = 1 ./ (1 + exp(k * (G - lower_threshold))) + 1 ./ (1 + exp(k * (upper_threshold - G)));
    
    A  = A  + dt * (D_A  * lap(A)  + alpha_A  * N ./ (1 + (R2 / K_R2A).^n) - delta_A  * A);
    G  = G  + dt * (D_G  * lap(G)  + alpha_G  * A - delta_G  * G);
    R  = R  + dt * (D_R  * lap(R)  + alpha_R  * N ./ (1 + (R / K_RR ).^n .* H_G) - delta_R  * R);
    R2 = R2 + dt * (D_R2 * lap(R2) + alpha_R2 * N ./ (1 + (R / K_RR2).^n .* H_G) - delta_R2 * R2);
    C  = C  + dt * (D_C  * lap(C)  + alpha_C  * N ./ (1 + (R2 / K_R2C).^n + (R / K_RC).^n .* H_G) - delta_C * C);
    
    % Dynamic Edge Tracking
    valid_pixels = N > N_edge_threshold;
    if any(valid_pixels(:))
        max_colony_radius = max(R_cells(valid_pixels));
        edge_mask = (R_cells <= max_colony_radius) & ...
                    (R_cells >= (max_colony_radius - ring_width)) & valid_pixels;
        
        N_edge(ti)  = mean(N(edge_mask));
        A_edge(ti)  = mean(A(edge_mask));
        G_edge(ti)  = mean(G(edge_mask));
        R_edge(ti)  = mean(R(edge_mask));
        R2_edge(ti) = mean(R2(edge_mask));
        C_edge(ti)  = mean(C(edge_mask));
    else
        N_edge(ti) = N(cy, cx); A_edge(ti) = A(cy, cx); G_edge(ti) = G(cy, cx);
        R_edge(ti) = R(cy, cx); R2_edge(ti) = R2(cy, cx); C_edge(ti) = C(cy, cx);
    end
end
toc

%  ==========================================
%  NORMALIZATION AND SUBPLOT GENERATION
%  ==========================================
% Normalize data to a maximum of 1
N_norm  = N_edge  / max(N_edge);
A_norm  = A_edge  / max(A_edge);
R_norm  = R_edge  / max(R_edge);
R2_norm = R2_edge / max(R2_edge);
C_norm  = C_edge  / max(C_edge);

%  ==========================================
%  CPK/N Normalization
%  ==========================================
% Ratios CPK/N at edge
C_N_edge = C_edge ./ N_edge;

% Create a tall figure for stacked subplots
fig = figure('Name', 'Stacked Oscillations at Colony Edge', 'Position', [100 100 800 1000], 'Color', 'w');

% Use tiledlayout to eliminate whitespace between plots and align axes perfectly
t = tiledlayout(6, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% 1. Population N
ax1 = nexttile;
plot(time_vector, N_edge, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]);
ylabel('Population N', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTickLabel', []); % Hide X-axis labels for all but the bottom plot
box on;

% 2. ScbA
ax2 = nexttile;
plot(time_vector, A_norm, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
ylabel('ScbA', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTickLabel', []);
box on;

% 3. ScbR
ax3 = nexttile;
plot(time_vector, R_norm, 'LineWidth', 2, 'Color', [0.0000 0.4470 0.7410]);
ylabel('ScbR', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTickLabel', []);
box on;

% 4. ScbR2
ax4 = nexttile;
plot(time_vector, R2_norm, 'LineWidth', 2, 'Color', [0.9290 0.6940 0.1250]);
ylabel('ScbR2', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTickLabel', []);
box on;

% 5. CPK
ax5 = nexttile;
plot(time_vector, C_norm, 'LineWidth', 2, 'Color', [0.4940 0.1840 0.5560]);
ylabel('CPK', 'FontSize', 12, 'FontWeight', 'bold');
set(gca, 'XTickLabel', []);
box on;

% 6. CPK/N
ax5 = nexttile;
plot(time_vector, C_N_edge, 'LineWidth', 3.5, 'Color', [0 0 0]);
ylabel('CPK/N', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Time', 'FontSize', 14, 'FontWeight', 'bold');
box on;

 %save CPK/N as a CSV file
 %CPK_N = [time_vector(:), C_N_edge(:)];
 %writematrix(CPK_N,'WT_model.csv')

% Link the x-axes so zooming/panning affects all plots simultaneously
linkaxes([ax1, ax2, ax3, ax4, ax5], 'x');
xlim(ax1, [0 max(time_vector)]);

% Set consistent Y-limits for all normalized plots
%set([ax1 ax2 ax3 ax4 ax5], 'YLim', [0 1.05], 'FontSize', 11, 'LineWidth', 1.2);

% Add an overarching title to the layout
title(t, sprintf('Simulated Network Oscillations at Colony Edge (Ring width = %.1f)', ring_width), 'FontSize', 16, 'FontWeight', 'bold');