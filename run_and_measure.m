% =========================================================
% run_and_measure.m
%
% Runs the calibrated QS-repressilator simulation once and measures
% the period and amplitude from the returned edge CPK/N trace.
%
% Requires in the same folder: run_qs_sim.m
% Requires: Signal Processing Toolbox (findpeaks)
% =========================================================
close all; clear; clc;

% ---------------------------------------------------------
% 1. INITIALIZATION
% ---------------------------------------------------------

OVERRIDES = struct();

% ---------------------------------------------------------
% 2. DETECTOR SETTINGS
% ---------------------------------------------------------

PROM_FRAC  = 0.20;    % prominence as a fraction of the trace range
MIN_SEP_H  = 8;       % minimum peak separation, hours
TRIM_H     = 0;       % discard this many hours of transient before measuring
                      % (0 = faithful to the engine and to the calibration)

% Amplitude pairing:
%   'legacy'  pairs the j-th peak with the j-th trough by index.
%             This is what run_qs_sim.m and the experimental script both do,
%             so it reproduces the 0.587 calibration convention.
%   'ordered' pairs each peak with the FIRST trough that follows it.
%             This is a true peak-to-valley drop.
% Both are printed regardless; this only sets which one is reported as primary.

PAIRING    = 'legacy';

SAVE_TRACE = false;   % write the measured trace to CSV
TRACE_FILE = 'sim_trace.csv';

% ---------------------------------------------------------
% 3. RUN
% ---------------------------------------------------------
fprintf('Running simulation...\n');
t0 = tic;
r  = run_qs_sim(OVERRIDES);
elapsed = toc(t0);

if r.blewup
    fprintf('\n*** Run went non-finite. No measurement possible. ***\n');
    return;
end

fprintf('Done in %.0f s.  Trace: %d samples over %.1f h.\n\n', ...
        elapsed, numel(r.t), r.t(end));

% ---------------------------------------------------------
% 4. MEASURE
% ---------------------------------------------------------
m = measure_trace(r.t, r.CN, PROM_FRAC, MIN_SEP_H, TRIM_H);

if m.n_peaks < 2
    fprintf('Only %d peak(s) detected. Not enough to define a period.\n', m.n_peaks);
    fprintf('Check the trace, the horizon, or lower PROM_FRAC.\n');
    return;
end

switch lower(PAIRING)
    case 'legacy',  amp_primary = m.amp_legacy;
    case 'ordered', amp_primary = m.amp_ordered;
    otherwise, error('PAIRING must be ''legacy'' or ''ordered''.');
end

% ---------------------------------------------------------
% 5. REPORT
% ---------------------------------------------------------
fprintf('=======================================================\n');
fprintf('MEASUREMENT\n');
fprintf('=======================================================\n');
fprintf('  peaks detected     : %d\n', m.n_peaks);
fprintf('  peak times (h)     : %s\n', mat2str(round(m.t_peaks',2)));
fprintf('  intervals (h)      : %s\n', mat2str(round(m.intervals',2)));
fprintf('\n');
fprintf('  PERIOD             : %.2f h\n', m.period);
if numel(m.intervals) > 1
    fprintf('    cycle-to-cycle sd: %.2f h  (within this one trace, %d intervals)\n', ...
            std(m.intervals,1), numel(m.intervals));
    fprintf('    NOTE: this is limit-cycle drift, NOT the experimental +/- 1.4 h.\n');
    fprintf('          A single deterministic run has no replicate spread.\n');
else
    fprintf('    only one interval, no drift estimate available\n');
end
fprintf('  period (autocorr)  : %.2f h  (engine cross-check)\n', r.period_acf);
fprintf('\n');
fprintf('  AMPLITUDE (%s) : %.3f\n', PAIRING, amp_primary);
fprintf('    legacy pairing   : %.3f  (%d cycles)\n', m.amp_legacy,  m.ncyc_legacy);
fprintf('    ordered pairing  : %.3f  (%d cycles)\n', m.amp_ordered, m.ncyc_ordered);
if abs(m.amp_legacy - m.amp_ordered) > 0.02
    fprintf('    *** The two pairings disagree by %.3f. The trace likely starts\n', ...
            abs(m.amp_legacy - m.amp_ordered));
    fprintf('        on a downswing, so the legacy index pairing is mixing\n');
    fprintf('        trough-to-peak rises in with peak-to-valley drops.\n');
end

% % ---------------------------------------------------------
% % 6. CROSS-CHECK AGAINST THE ENGINE
% % ---------------------------------------------------------
% 
% fprintf('\n-------------------------------------------------------\n');
% fprintf('ENGINE CROSS-CHECK (run_qs_sim internal measurement)\n');
% fprintf('  engine period      : %.2f h   here: %.2f h   diff %.3f\n', ...
%         r.period_pk, m.period, abs(r.period_pk - m.period));
% fprintf('  engine amplitude   : %.3f     here: %.3f     diff %.3f\n', ...
%         r.amp, m.amp_legacy, abs(r.amp - m.amp_legacy));
% if TRIM_H == 0
%     if abs(r.period_pk - m.period) < 1e-6 && abs(r.amp - m.amp_legacy) < 1e-6
%         fprintf('  MATCH. This script reproduces the engine exactly.\n');
%     else
%         fprintf('  MISMATCH with TRIM_H = 0. Detector settings here differ\n');
%         fprintf('  from the ones hard-coded in run_qs_sim.m lines 122-131.\n');
%     end
% else
%     fprintf('  (TRIM_H = %.1f h, so a difference from the engine is expected.)\n', TRIM_H);
% end
% fprintf('=======================================================\n');

% ---------------------------------------------------------
% 7. PLOT
% ---------------------------------------------------------
figure('Color','w','Position',[100 100 950 420]); hold on; box on;
plot(m.t, m.sig_n, '-', 'Color',[0 0.447 0.741], 'LineWidth',1.5);
plot(m.t(m.pk), m.sig_n(m.pk), 'v', 'Color',[0.85 0 0], ...
     'MarkerFaceColor',[0.85 0 0], 'MarkerSize',7);
plot(m.t(m.tr), m.sig_n(m.tr), '^', 'Color',[0 0.5 0], ...
     'MarkerFaceColor',[0 0.5 0], 'MarkerSize',7);
xlabel('Time (h)'); ylabel('Normalized CPK/N  [0,1]');
title(sprintf('period %.2f h   |   amplitude %.3f (%s)   |   %d peaks', ...
      m.period, amp_primary, PAIRING, m.n_peaks));
xlim([m.t(1) m.t(end)]); ylim([-0.05 1.05]);
legend({'CPK/N (normalized)','peaks','troughs'}, 'Location','best');
hold off;

% ---------------------------------------------------------
% 8. OPTIONAL EXPORT
% ---------------------------------------------------------
if SAVE_TRACE
    writematrix([m.t, m.sig_raw, m.sig_n], TRACE_FILE);
    fprintf('\nTrace written to %s  [t_h, CPK_N_raw, CPK_N_norm]\n', TRACE_FILE);
end


% =========================================================
% LOCAL FUNCTION
% =========================================================
function m = measure_trace(tt, sig, prom_frac, min_sep_h, trim_h)
% Measure period and amplitude from a single time series.

    keep = tt >= (tt(1) + trim_h);
    tt   = tt(keep);
    sig  = sig(keep);

    dt_h  = tt(2) - tt(1);
    mind  = round(min_sep_h / dt_h);

    sig_n = (sig - min(sig)) / (max(sig) - min(sig));
    rng_n = max(sig_n) - min(sig_n);

    [~, pk] = findpeaks( sig_n, 'MinPeakProminence', prom_frac*rng_n, ...
                                'MinPeakDistance',   mind);
    [~, tr] = findpeaks(-sig_n, 'MinPeakProminence', prom_frac*rng_n, ...
                                'MinPeakDistance',   mind);

    % ---- period ----
    if numel(pk) >= 2
        intervals = diff(tt(pk));
        period    = mean(intervals);
    else
        intervals = [];
        period    = NaN;
    end

    % ---- amplitude, legacy index pairing (matches run_qs_sim.m) ----
    ncyc_l = min(numel(pk), numel(tr));
    if ncyc_l >= 1
        amp_l = mean(sig_n(pk(1:ncyc_l)) - sig_n(tr(1:ncyc_l)));
    else
        amp_l = NaN;
    end

    % ---- amplitude, ordered pairing (peak to the next trough after it) ----
    drops = [];
    for a = 1:numel(pk)
        nxt = tr(tr > pk(a));
        if ~isempty(nxt)
            drops(end+1,1) = sig_n(pk(a)) - sig_n(nxt(1)); %#ok<AGROW>
        end
    end
    if isempty(drops), amp_o = NaN; else, amp_o = mean(drops); end

    m = struct( ...
        't',            tt, ...
        'sig_raw',      sig, ...
        'sig_n',        sig_n, ...
        'pk',           pk, ...
        'tr',           tr, ...
        't_peaks',      tt(pk), ...
        'intervals',    intervals, ...
        'period',       period, ...
        'n_peaks',      numel(pk), ...
        'amp_legacy',   amp_l, ...
        'ncyc_legacy',  ncyc_l, ...
        'amp_ordered',  amp_o, ...
        'ncyc_ordered', numel(drops));
end
