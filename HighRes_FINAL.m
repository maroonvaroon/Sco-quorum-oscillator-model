close all; clear;

tic

%  =======================
%  SPEED / RESOLUTION
%  =======================

dx_base = 1;
dt_base = 0.1;

dx_ref  = 0.125;
T_ref = 49000;
                                                                                                  
tau_cal = 0.70597;

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
Nx = max(4, round(Lx / dx));
Ny = max(4, round(Ly / dx));


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
eps = 0.01; % Sharpness of the border (smaller=larger)

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
%  COLOR CONTROL
%  ==========================

% Data kymograph color palette
colors = [0 0 0        % black
    0.6275 0.1255 0.9412 % purple
    0 0 1                % blue
    0 1 1                % cyan
    0 1 0                % green
    1 1 0                % yellow
    1 0.6471 0           % orange
    1 0 0                % red
    1 0.7529 0.7961]     % pink

% interpolate to smooth the color scheme
ncolors = 100
x0 = linspace(0,1,size(colors,1));
xf = linspace(0,1,ncolors);
CustomPalette = interp1(x0,colors,xf)

%  ==========================
%  INITIAL CONDITIONS
%  ==========================

% Initial Position

[X, Y] = meshgrid(1:Nx, 1:Ny);
x0 = Nx / 2;
y0 = Ny / 2;

% Initial population spread/dipersion (Gaussian)

sigma_phys = 3; % (Larger -> More dispersed)

sigma = sigma_phys / dx;

N  = K_N * 0.01 * exp(-((X - x0).^2 + (Y - y0).^2) / (2 * sigma^2));
A  = zeros(Ny, Nx);
G  = A;
R  = A;
R2 = A;
C  = A;
A(round(y0), round(x0)) = 0.1;

%  =================================
%  Laplacian operator for the Reaction/Diffusion model
%  =================================

if USE_LAPLACIAN_9
    lap = @(U) ( ...
        0.25 * (circshift(U,[1,1]) + circshift(U,[-1,-1]) + ...
                circshift(U,[1,-1]) + circshift(U,[-1,1])) + ...
        0.5  * (circshift(U,[1,0]) + circshift(U,[-1,0]) + ...
                circshift(U,[0,1]) + circshift(U,[0,-1])) - ...
        3 * U ) / dx^2;
else
    lap = @(U) ( ...
        circshift(U,[1,0]) + circshift(U,[-1,0]) + ...
        circshift(U,[0,1]) + circshift(U,[0,-1]) - ...
        4 * U ) / dx^2;
end

%  ==========================================
%  KYMOGRAPH (radial bins)
%  ==========================================
R_cells = hypot((X - x0), (Y - y0));      

R_MAX_LEN = 50;                          
r_max_full = floor(min(Nx, Ny)/2);       
r_max = min(r_max_full, floor(R_MAX_LEN / dx)); 

r_centers = 0:r_max;                      
edges     = (-0.5:1:(r_max+0.5));         

bin       = discretize(R_cells, edges);   
Nr        = numel(r_centers);

validMask = ~isnan(bin);
cntPerBin = accumarray(bin(validMask), 1, [Nr 1], @sum, 0);

nKymoRows  = ceil(T / KYMO_EVERY);
kymo_C     = nan(nKymoRows, Nr);
kymo_N     = nan(nKymoRows, Nr);
kymo_CNpx  = nan(nKymoRows, Nr);
kymo_t     = nan(nKymoRows, 1);

r_cent_len = r_centers * dx;             

%  =========================
%  FIGURES
%  =========================
fig = figure('Name','Concentrations + Kymographs','Position',[60 60 1600 950]);

tlo = tiledlayout(fig, 5, 3, 'TileSpacing','compact', 'Padding','compact');

axN  = nexttile(tlo, 1); hN  = imagesc(N,  'Parent', axN); axis(axN,'square','off'); colormap(axN,CustomPalette); colorbar(axN); title(axN,'Population N');
axA  = nexttile(tlo, 2); hA  = imagesc(A,  'Parent', axA); axis(axA,'square','off'); colormap(axA,CustomPalette); colorbar(axA); title(axA,'ScbA');
axG  = nexttile(tlo, 3); hG  = imagesc(G,  'Parent', axG); axis(axG,'square','off'); colormap(axG,CustomPalette); colorbar(axG); title(axG,'GBL');
axR  = nexttile(tlo, 4); hR  = imagesc(R,  'Parent', axR); axis(axR,'square','off'); colormap(axR,CustomPalette); colorbar(axR); title(axR,'ScbR');
axR2 = nexttile(tlo, 5); hR2 = imagesc(R2, 'Parent', axR2);axis(axR2,'square','off'); colormap(axR2,CustomPalette); colorbar(axR2); title(axR2,'ScbR2');
axC  = nexttile(tlo, 6); hC  = imagesc(C,  'Parent', axC); axis(axC,'square','off'); colormap(axC,CustomPalette); colorbar(axC); title(axC,'CPK');

axK = nexttile(tlo, 7,  [1 3]);
hK  = imagesc(r_cent_len, 0, nan(1, Nr), 'Parent', axK);
set(axK,'YDir','reverse'); colormap(axK,CustomPalette); colorbar(axK);
xlabel(axK,'Radius (length units)'); ylabel(axK,'Time'); title(axK,'CPK Radial Kymograph');
caxis(axK,[0 7]); set(axK,'YLim',[0 eps]);

axKN = nexttile(tlo, 10, [1 3]);
hKN  = imagesc(r_cent_len, 0, nan(1, Nr), 'Parent', axKN);
set(axKN,'YDir','reverse'); colormap(axKN,CustomPalette); colorbar(axKN);
xlabel(axKN,'Radius (length units)'); ylabel(axKN,'Time'); title(axKN,'Population N Radial Kymograph');
caxis(axKN,[0 1.1*K_N]); set(axKN,'YLim',[0 eps]);

axKCN = nexttile(tlo, 13, [1 3]);
hKCN  = imagesc(r_cent_len, 0, nan(1, Nr), 'Parent', axKCN);
set(axKCN,'YDir','reverse'); colormap(axKCN,CustomPalette); colorbar(axKCN);
xlabel(axKCN,'Radius (length units)'); ylabel(axKCN,'Time'); title(axKCN,'CPK/N Radial Kymograph');
caxis(axKCN, [0 5.5]); set(axKCN,'YLim',[0 eps]);

xlim(axK,   [0 R_MAX_LEN]);
xlim(axKN,  [0 R_MAX_LEN]);
xlim(axKCN, [0 R_MAX_LEN]);

drawnow;

%  ==========================
%  POINTWISE CPK/N FIELD
%  ==========================
N_min_px = 1e-3 * K_N;          % same occupancy threshold as the kymograph

figCN = figure('Name','Pointwise CPK/N','Position',[1700 60 620 560],'Color','k');
axCN  = axes(figCN);
hCN   = imagesc(nan(Ny,Nx), 'Parent', axCN);
set(hCN, 'AlphaData', zeros(Ny,Nx));   % un-colonised pixels stay blank
axis(axCN,'square','off');
colormap(axCN, CustomPalette); colorbar(axCN);

clim(axCN, [0 5.5]);     %fixes color scale

title(axCN,'CPK/N');

drawnow;

%% ================
%  TIME LOOP
%  ================
kymoRow = 0;

N_min  = 1e-3 * K_N;
N_safe = N_min;

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

    % Bandpass Filter (k is the sharpness of the transition)

    H_G = 1 ./ (1 + exp(k * (G - lower_threshold))) + ...
          1 ./ (1 + exp(k * (upper_threshold - G)));

    A  = A  + dt * (D_A  * lap(A)  + alpha_A  * N ./ (1 + (R2 / K_R2A).^n) - delta_A  * A);
    G  = G  + dt * (D_G  * lap(G)  + alpha_G  * A - delta_G  * G);
    R  = R  + dt * (D_R  * lap(R)  + alpha_R  * N ./ (1 + (R / K_RR ).^n .* H_G) - delta_R  * R);
    R2 = R2 + dt * (D_R2 * lap(R2) + alpha_R2 * N ./ (1 + (R / K_RR2).^n .* H_G) - delta_R2 * R2);
    C  = C  + dt * (D_C  * lap(C)  + alpha_C  * N ./ ...
          (1 + (R2 / K_R2C).^n + (R / K_RC).^n .* H_G) - delta_C * C);

    if mod(ti, KYMO_EVERY) == 0
        kymoRow = kymoRow + 1;

        sumsC = accumarray(bin(validMask), C(validMask), [Nr 1], @sum, 0);
        sumsN = accumarray(bin(validMask), N(validMask), [Nr 1], @sum, 0);

        kymo_C(kymoRow,:) = (sumsC ./ max(cntPerBin,1)).';
        kymo_N(kymoRow,:) = (sumsN ./ max(cntPerBin,1)).';

        occ = (N > N_min);
        CNpx = zeros(size(C));
        CNpx(occ) = C(occ) ./ max(N(occ), N_safe);

        sumsCN = accumarray(bin(validMask), CNpx(validMask), [Nr 1], @sum, 0);
        kymo_CNpx(kymoRow,:) = (sumsCN ./ max(cntPerBin,1)).';

        kymo_t(kymoRow) = (ti-1)*dt*tau_cal;
    end

    if mod(ti, PLOT_EVERY) == 0
        set(hN,  'CData', N);
        set(hA,  'CData', A);
        set(hG,  'CData', G);
        set(hR,  'CData', R);
        set(hR2, 'CData', R2);
        set(hC,  'CData', C);
    end

    if mod(ti, KYMO_EVERY) == 0
        set(hK,   'CData', kymo_C(1:kymoRow, :),    'XData', r_cent_len, 'YData', kymo_t(1:kymoRow));
        set(hKN,  'CData', kymo_N(1:kymoRow, :),    'XData', r_cent_len, 'YData', kymo_t(1:kymoRow));
        set(hKCN, 'CData', kymo_CNpx(1:kymoRow, :), 'XData', r_cent_len, 'YData', kymo_t(1:kymoRow));

        yMax = max(kymo_t(kymoRow), eps);
        set(axK,   'YLim', [0, yMax]);
        set(axKN,  'YLim', [0, yMax]);
        set(axKCN, 'YLim', [0, yMax]);
    end

    % ---- pointwise CPK/N ----
    if mod(ti, PLOT_EVERY) == 0
        occ_px = N > N_min_px;
        CN_px  = nan(Ny, Nx);
        CN_px(occ_px) = C(occ_px) ./ N(occ_px);

        set(hCN, 'CData', CN_px, 'AlphaData', occ_px);

        hi = max(CN_px(occ_px));
        if isempty(hi) || ~isfinite(hi) || hi <= 0, hi = 1; end
        caxis(axCN, [0 5.5]);

        title(axCN, sprintf('CPK/N (pointwise), t = %.2f', (ti-1)*dt*tau_cal));
    end

    if mod(ti, PLOT_EVERY) == 0 || mod(ti, KYMO_EVERY) == 0
        drawnow limitrate;
    end
end

% Final refresh
set(hK,   'CData', kymo_C(1:kymoRow, :),    'XData', r_cent_len, 'YData', kymo_t(1:kymoRow));
set(hKN,  'CData', kymo_N(1:kymoRow, :),    'XData', r_cent_len, 'YData', kymo_t(1:kymoRow));
set(hKCN, 'CData', kymo_CNpx(1:kymoRow, :), 'XData', r_cent_len, 'YData', kymo_t(1:kymoRow));

yMax = max(kymo_t(kymoRow), eps);
set(axK,   'YLim', [0, yMax]);
set(axKN,  'YLim', [0, yMax]);
set(axKCN, 'YLim', [0, yMax]);
drawnow;

toc
