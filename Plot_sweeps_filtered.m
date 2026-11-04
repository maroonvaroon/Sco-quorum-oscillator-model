%% Plot_sweeps_filtered.m

clear; close all; clc;

%% ---------------------------------------------------------------- options
LAYOUT = 'all';        % 'DG_KN' (2x2) | 'all' (2x7) | 'all_portrait' (7x2)
STYLE  = 'errorbar';   % 'errorbar' | 'shaded'

% Point Highlight

HIGHLIGHT.DG  = 3.0;
HIGHLIGHT.K_N = 1000;

% HIGHLIGHT.mu      = [];
% HIGHLIGHT.alpha_G = [];
% HIGHLIGHT.delta_G = [];
% HIGHLIGHT.K_RC    = [];
% HIGHLIGHT.K_R2C   = [];

OUT_STEM = ['sweeps_filtered_' LAYOUT '_' STYLE];

%% ----------------------------------------------------------- Filter Parameters
TC         = 6;      % h, cutoff period
ORDER      = 3;      % Butterworth order
MINSEP     = 5;      % h, minimum peak/trough separation
PROM       = 0.15;   % prominence as a fraction of the filtered-trace range
T_START    = 9;      % h, extrema before this are ignored
MIN_CYCLES = 2;      % runs with fewer intervals or cycles are omitted

PERIOD_YLIM = [7 20];
AMP_YLIM    = [0 1];

%% ------------------------------------------------------- Time Series
% name, file, variable holding the sweep values, axis label, nominal value,
% x-scale factor for display, x-label, [vmin vmax] range to keep
sw = {};
sw(end+1,:) = {'DG',      'DG_sweep_results.mat',      'DG_VALUES', 'D_G',          1.0,  1,    'D_G',                [-Inf Inf]};
sw(end+1,:) = {'mu',      'mu_sweep_results.mat',      'VALUES',    '\mu',          0.3,  1,    '\mu',                [0.2 0.5]};
sw(end+1,:) = {'alpha_G', 'alpha_G_sweep_results.mat', 'VALUES',    '\alpha_G',     1.0,  1,    '\alpha_G',           [-Inf Inf]};
sw(end+1,:) = {'delta_G', 'delta_G_sweep_results.mat', 'VALUES',    '\delta_G',     1.0,  1,    '\delta_G',           [-Inf Inf]};
sw(end+1,:) = {'K_RC',    'K_RC_sweep_results.mat',    'VALUES',    'K_{RC}',       1.0,  1,    'K_{RC}',             [-Inf Inf]};
sw(end+1,:) = {'K_R2C',   'K_R2C_sweep_results.mat',   'VALUES',    'K_{R2C}',      1.0,  1,    'K_{R2C}',            [-Inf Inf]};
sw(end+1,:) = {'K_N',     'K_N_sweep_results.mat',     'VALUES',    'K_N',          2000, 1e-3, 'K_N (\times 10^3)',  [-Inf Inf]};

sweeps = cell2struct(sw, {'name','file','vkey','label','nominal','xscale','xlabel','vrange'}, 2);

if strcmp(LAYOUT, 'DG_KN')
    sweeps = sweeps(ismember({sweeps.name}, {'DG', 'K_N'}));
end
nS = numel(sweeps);

%% ------------------------------------------------------------ measurement
for s = 1:nS
    d = load(sweeps(s).file);
    x = d.(sweeps(s).vkey)(:);
    T = measure_sweep(d.traces, TC, ORDER, MINSEP, PROM, T_START);
    T.value = x;
    T = T(:, [end, 1:end-1]);

    keep = T.n_intervals >= MIN_CYCLES & T.n_cycles >= MIN_CYCLES ...
         & T.value >= sweeps(s).vrange(1) & T.value <= sweeps(s).vrange(2);
    T = T(keep, :);

    sweeps(s).T = T;
    writetable(T, [sweeps(s).name '_filtered_period_amplitude.csv']);

    iw = find(abs(T.value - sweeps(s).nominal) < 1e-9, 1);
    if ~isempty(iw)
        fprintf('%-8s WT: period %.2f +/- %.2f h, amplitude %.3f +/- %.3f\n', ...
            sweeps(s).name, T.period_mean(iw), T.period_sd(iw), T.amp_mean(iw), T.amp_sd(iw));
    end
end

%% ------------------------------------------------------------------ figure
switch LAYOUT
    case 'DG_KN'
        fig = figure('Units', 'inches', 'Position', [1 1 7.2 5.2], 'Color', 'w');
        tl  = tiledlayout(2, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
        tileOf = @(r, c) (r - 1) * nS + c;        % rows = metric, cols = parameter
    case 'all'
        fig = figure('Units', 'inches', 'Position', [1 1 16 5.2], 'Color', 'w');
        tl  = tiledlayout(2, nS, 'TileSpacing', 'compact', 'Padding', 'compact');
        tileOf = @(r, c) (r - 1) * nS + c;
    case 'all_portrait'
        fig = figure('Units', 'inches', 'Position', [1 1 7.2 15], 'Color', 'w');
        tl  = tiledlayout(nS, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
        tileOf = @(r, c) (c - 1) * 2 + r;         % rows = parameter, cols = metric
    otherwise
        error('LAYOUT must be ''DG_KN'', ''all'' or ''all_portrait''.');
end
portrait = strcmp(LAYOUT, 'all_portrait');

metrics = {'period', 'Period (h)', PERIOD_YLIM; 'amp', 'Amplitude', AMP_YLIM};

for s = 1:nS
    T   = sweeps(s).T;
    x   = T.value * sweeps(s).xscale;
    isNom = abs(T.value - sweeps(s).nominal) < 1e-9;
    if isfield(HIGHLIGHT, sweeps(s).name) && ~isempty(HIGHLIGHT.(sweeps(s).name))
        isHi = ismembertol(T.value, HIGHLIGHT.(sweeps(s).name), 1e-9, 'DataScale', 1);
    else
        isHi = false(size(x));
    end

    for r = 1:2
        ax = nexttile(tl, tileOf(r, s));
        mu = T.([metrics{r,1} '_mean']);
        sd = T.([metrics{r,1} '_sd']);
        draw_panel(ax, x, mu, sd, isNom, isHi, STYLE);
        ylim(ax, metrics{r,3});

        if portrait
            xlabel(ax, sweeps(s).xlabel);
            if s == 1, title(ax, metrics{r,2}); end
            if r == 1, ylabel(ax, sweeps(s).label, 'Rotation', 0, 'HorizontalAlignment', 'right'); end
        else
            if r == 1, title(ax, sweeps(s).label); end
            if r == 2, xlabel(ax, sweeps(s).xlabel); end
            if s == 1
                ylabel(ax, metrics{r,2});
            else
                ax.YTickLabel = [];
            end
        end
    end
end

% legends
hold(ax, 'on');
hL = gobjects(0, 1);
switch STYLE
    case 'errorbar'
        hL(end+1) = errorbar(ax, NaN, NaN, NaN, 'o-', 'Color', [0 0.447 0.741], ...
            'MarkerFaceColor', [0 0.447 0.741], 'MarkerSize', 3, 'DisplayName', 'mean \pm 1 SD');
    case 'shaded'
        hL(end+1) = fill(ax, NaN(1,3), NaN(1,3), [0 0.447 0.741], 'FaceAlpha', 0.2, ...
            'EdgeColor', 'none', 'DisplayName', '\pm 1 SD');
end
hL(end+1) = plot(ax, NaN, NaN, 'o', 'Color', [0.85 0.15 0.15], 'MarkerFaceColor', [0.85 0.15 0.15], ...
    'MarkerSize', 6, 'DisplayName', 'WT (nominal)');
hL(end+1) = plot(ax, NaN, NaN, 'o', 'Color', [0 0.5 0], 'MarkerFaceColor', [0 0.5 0], ...
    'MarkerSize', 6, 'DisplayName', 'highlighted');
lg = legend(hL, 'Orientation', 'horizontal', 'Box', 'off');
lg.Layout.Tile = 'north';

exportgraphics(fig, [OUT_STEM '.pdf'], 'ContentType', 'vector');
exportgraphics(fig, [OUT_STEM '.png'], 'Resolution', 300);
fprintf('Saved %s.pdf / .png\n', OUT_STEM);

%% ================================================================ functions
function T = measure_sweep(traces, Tc, order, minsep, prom, tStart)
    n = numel(traces);
    [pm, psd, am, asd] = deal(nan(n, 1));
    [ni, nc]           = deal(zeros(n, 1));
    for k = 1:n
        t  = traces{k}(:, 1);
        y  = traces{k}(:, 2);
        yf = lowpass_padded(t, y, Tc, order);
        [pk, tr] = extrema(t, yf, minsep, prom, tStart);

        iv = [diff(t(pk)); diff(t(tr))];

        sw = zeros(0, 1);
        for j = pk'
            nxt = tr(find(tr > j, 1));
            if ~isempty(nxt)
                sw(end + 1, 1) = (yf(j) - yf(nxt)) / max(yf); %#ok<AGROW>
            end
        end

        ni(k) = numel(iv);  nc(k) = numel(sw);
        if ni(k) > 0, pm(k) = mean(iv); end
        if ni(k) > 1, psd(k) = std(iv); end
        if nc(k) > 0, am(k) = mean(sw); end
        if nc(k) > 1, asd(k) = std(sw); end
    end
    T = table(pm, psd, ni, am, asd, nc, 'VariableNames', ...
        {'period_mean', 'period_sd', 'n_intervals', 'amp_mean', 'amp_sd', 'n_cycles'});
end

function yf = lowpass_padded(t, y, Tc, order)
    dt    = mean(diff(t));
    [b,a] = butter(order, (1 / Tc) / (0.5 / dt), 'low');
    y     = y(:);
    np    = min(floor(3 * Tc / dt), numel(y) - 1);
    ye    = [2*y(1) - y(np+1:-1:2); y; 2*y(end) - y(end-1:-1:end-np)];
    yf    = filtfilt(b, a, ye);
    yf    = yf(np + 1 : np + numel(y));
end

function [pk, tr] = extrema(t, y, minsep, prom, tStart)
    dt = mean(diff(t));
    p  = prom * (max(y) - min(y));
    d  = floor(minsep / dt);
    [~, pk] = findpeaks( y, 'MinPeakProminence', p, 'MinPeakDistance', d);
    [~, tr] = findpeaks(-y, 'MinPeakProminence', p, 'MinPeakDistance', d);
    pk = pk(:);  pk = pk(t(pk) >= tStart);
    tr = tr(:);  tr = tr(t(tr) >= tStart);
end

function draw_panel(ax, x, mu, sd, isNom, isHi, style)
    blue  = [0 0.447 0.741];
    red   = [0.85 0.15 0.15];
    green = [0.2 1 0.2];
    x = x(:); mu = mu(:); sd = sd(:);
    hold(ax, 'on');
    switch style
        case 'errorbar'
            errorbar(ax, x, mu, sd, 'o-', 'Color', blue, 'MarkerFaceColor', blue, ...
                'MarkerSize', 3, 'LineWidth', 0.8, 'CapSize', 2);
            if any(isHi)
                errorbar(ax, x(isHi), mu(isHi), sd(isHi), 'o', 'Color', green, ...
                    'MarkerFaceColor', green, 'MarkerSize', 6, 'LineWidth', 1.2, 'CapSize', 3);
            end
            w = isNom & ~isHi;
            if any(w)
                errorbar(ax, x(w), mu(w), sd(w), 'o', 'Color', 'magenta', ...
                    'MarkerFaceColor', red, 'MarkerSize', 6, 'LineWidth', 1.2, 'CapSize', 3);
            end
        case 'shaded'
            ok = ~isnan(sd);
            fill(ax, [x(ok); flipud(x(ok))], [mu(ok) - sd(ok); flipud(mu(ok) + sd(ok))], blue, ...
                'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(ax, x, mu, '-', 'Color', blue, 'LineWidth', 1.8);
            if any(isHi)
                plot(ax, x(isHi), mu(isHi), 'o', 'Color', green, 'MarkerFaceColor', green, 'MarkerSize', 6.5);
            end
            w = isNom & ~isHi;
            if any(w)
                plot(ax, x(w), mu(w), 'o', 'Color', red, 'MarkerFaceColor', red, 'MarkerSize', 6.5);
            end
        otherwise
            error('STYLE must be ''errorbar'' or ''shaded''.');
    end
    pad = 0.03 * (max(x) - min(x));
    xlim(ax, [min(x) - pad, max(x) + pad]);
    grid(ax, 'on');
    box(ax, 'on');
end
