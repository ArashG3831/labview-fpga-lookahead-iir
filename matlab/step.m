% ===================== Step-response check (timestamp-only) =====================
% Compare FPGA filtered step vs MATLAB IIR with ONLY:
%   - timestamp alignment
%   - known fixed streaming delay of -4 samples (for this host link)
% No cross-correlation, no sign flip, no gain/offset fit.

clear; clc; close all;

% ----------- Inputs (two columns: time, value) -----------
raw_file  = 'raw4.txt';   % pre-filter ADC stream (step input): [t, x]
fpga_file = 'fi4.txt';    % FPGA-filtered stream: [t, y_fpga]

% Your look-ahead biquad (k=1)
b = [0.04613180209331293, 0.15257101841705373,...
     0.18943481534878998, 0.10568378381967050, 0.02268818479462131];
a = [1, 0, -0.7253696722084265, 0, 0.24187927668187462];

% Known deterministic pipeline delay on THIS step capture:
PIPE_DELAY = -4;  % samples (negative ⇒ FPGA output lags MATLAB by 4 samples)

% ---------------------- Load & stitch (timestamp only) -------------------
[t,  x_raw,  Fs]  = load_trace_stitched(raw_file);
[tf, y_hw,   FsF] = load_trace_stitched(fpga_file);
assert(abs(Fs-FsF) <= 1e-9*Fs, 'Fs mismatch: raw=%.9g, fpga=%.9g', Fs, FsF);
N  = numel(t);
fprintf('Fs = %.6f MHz, N=%d  (PIPE_DELAY = %d samples)\n', Fs/1e6, N, PIPE_DELAY);

% Align streams by integer offset from their start timestamps (no xcorr)
start_offset_samples = round((tf(1) - t(1)) * Fs);
y_on_raw_grid        = circshift(y_hw, -start_offset_samples);

% Apply ONLY the known fixed streaming delay
y_fpga = circshift(y_on_raw_grid, PIPE_DELAY);

% ---------------------- MATLAB reference (zero ICs) ----------------------
y_ideal = filter(b, a, double(x_raw));

% ---------------------- Find the step and compute metrics ----------------
% Robust step index from derivative + MAD threshold
dx   = [0; diff(double(x_raw))];
madx = median(abs(dx - median(dx))) + eps;
thr  = 6*madx;                                  % conservative
cand = find(abs(dx) > thr, 1, 'first');
if isempty(cand), error('No step edge detected in raw input.'); end
step_idx = cand;         % first significant edge
t_step   = t(step_idx);

% Define pre/post windows for statistics
preN  = min(step_idx-5, max(100, round(0.02*N)));
postN = min(N-step_idx-5, max(500,  round(0.05*N)));
preIdx  = (step_idx-preN):(step_idx-5);
postIdx = (step_idx+5):(step_idx+postN);

x0 = median(double(x_raw(preIdx)));
x1 = median(double(x_raw(postIdx)));
u_step = x1 - x0;                               % input step amplitude

% Output baselines/finals
yi0 = median(y_ideal(preIdx));    yi1 = median(y_ideal(postIdx));
yf0 = median(y_fpga (preIdx));    yf1 = median(y_fpga (postIdx));

% DC gain (measured) and theoretical from coefficients
G_ideal = (yi1 - yi0) / max(u_step,eps);
G_fpga  = (yf1 - yf0) / max(u_step,eps);
G_theory = sum(b) / sum(a);      % filter() convention ⇒ DC gain = sum(b)/sum(a)

% Step metrics (10–90% rise, overshoot, 1% settling)
mI = step_metrics(y_ideal, t, step_idx, yi0, yi1);
mF = step_metrics(y_fpga,  t, step_idx, yf0, yf1);

% ---------------------- Print summary ------------------------------------
fprintf('\n=== Step metrics (timestamp-only, no post-fit) ===\n');
fprintf('Input step amplitude Δu = %.6g\n', u_step);
fprintf('Measured DC gain:  Ideal=%.6g,  FPGA=%.6g,  Theory=%.6g\n', G_ideal, G_fpga, G_theory);

disp('Ideal (MATLAB):');
fprintf('  Rise time (10–90%%): %.3f ns\n', 1e9*mI.rise_10_90);
fprintf('  Peak overshoot: %.3f %% at t = %.3f ns\n', 100*mI.overshoot_frac, 1e9*mI.t_peak);
fprintf('  Settling time (±1%%): %.3f ns\n', 1e9*mI.settle_1pct);

disp('FPGA:');
fprintf('  Rise time (10–90%%): %.3f ns\n', 1e9*mF.rise_10_90);
fprintf('  Peak overshoot: %.3f %% at t = %.3f ns\n', 100*mF.overshoot_frac, 1e9*mF.t_peak);
fprintf('  Settling time (±1%%): %.3f ns\n', 1e9*mF.settle_1pct);

% Error over a steady window after the step
win0 = step_idx + round(0.02*N);                % start a bit after the edge
win1 = min(N - round(0.01*N), step_idx + postN);
win  = win0:win1;
e        = y_fpga(win) - y_ideal(win);
rel_rms  = rms(e) / max(rms(y_ideal(win)), eps);
nmse_db  = 20*log10(rel_rms);
fprintf('Relative RMS error (post-step) = %.4g  (NMSE = %.2f dB)\n\n', rel_rms, nmse_db);

% ---------------------- Plots --------------------------------------------
tt = (t - t(1))*1e6;

% 1) Overlaid step responses
figure('Color','w','Name','Step response — FPGA vs MATLAB (timestamp + fixed delay)');
plot(tt, y_ideal, 'k','LineWidth',1.1); hold on;
plot(tt, y_fpga , 'r--','LineWidth',1.0);
xline((t_step - t(1))*1e6, ':', 'Step', 'LabelVerticalAlignment','bottom');
grid on; xlabel('Time (\mus)'); ylabel('Amplitude (a.u.)');
title(sprintf('Filters only — FPGA (unaltered) vs MATLAB  [delay=%+d]', PIPE_DELAY));
legend('MATLAB ideal','FPGA (timestamp + fixed delay)','Location','best');

% 2) Zoom around the edge
zoomSpan = max(200, round(0.05*N));        % samples around the edge
z0 = max(1, step_idx - zoomSpan);
z1 = min(N, step_idx + zoomSpan);
figure('Color','w','Name','Step edge (zoom)');
plot(tt(z0:z1), y_ideal(z0:z1), 'k','LineWidth',1.1); hold on;
plot(tt(z0:z1), y_fpga (z0:z1), 'r--','LineWidth',1.0);
xline((t_step - t(1))*1e6, ':', 'Step');
grid on; xlabel('Time (\mus)'); ylabel('Amplitude (a.u.)');
title('Step edge (zoom)');

% 3) Error trace (post-step)
figure('Color','w','Name','Post-step error (time)');
plot(tt(win), e, 'LineWidth',1.0); grid on;
xlabel('Time (\mus)'); ylabel('Error (FPGA - MATLAB)');
title('Error after step (no post-fit)');

% ============================== Helpers =================================
function [t, x, Fs] = load_trace_stitched(fname)
    if ~exist(fname,'file'), error('File not found: %s', fname); end
    try
        T = readmatrix(fname); A = T(:,1:2);
    catch
        fid=fopen(fname,'r'); C=textscan(fid,'%f%f','CollectOutput',true); fclose(fid);
        A = C{1};
    end
    A = A(all(isfinite(A),2),:);
    t0 = A(:,1); x = A(:,2);

    % infer dt from positive diffs
    dt = diff(t0); pos = dt(dt>0);
    meddt = median(pos);
    if isempty(meddt) || ~isfinite(meddt), error('Cannot infer dt from %s', fname); end

    % stitch segments (time resets)
    t = t0; offset = 0; prev = t0(1);
    for i=2:numel(t0)
        if t0(i) <= prev, offset = offset + prev + meddt; end
        t(i) = t0(i) + offset; prev = t0(i);
    end
    Fs = 1/meddt;
end

function M = step_metrics(y, t, k0, y0, y1)
    % Compute basic step metrics for response y given step at index k0.
    y = double(y(:)); t = t(:);
    A  = y1 - y0;                     % expected change
    sgn = sign(A) + (A==0);           % handle zero safely
    y10 = y0 + 0.10*A;  y90 = y0 + 0.90*A;

    % First crossings after the step
    k10 = first_crossing(y, k0, y10, sgn);
    k90 = first_crossing(y, k10, y90, sgn);
    t10 = interp_time(t, y, k10, y10);
    t90 = interp_time(t, y, k90, y90);

    % Peak overshoot (signed)
    [pk, kp] = max(sgn*(y(k0:end)-y1));  % search after step
    kp = kp + k0 - 1;
    overshoot_frac = sgn*(y(kp) - y1) / max(abs(A), eps);
    t_peak = t(kp);

    % Settling time to ±1%
    band = 0.01*abs(A) + eps;
    ks = k90;   % start checking after 90% point
    settle_idx = NaN;
    for k = ks:numel(y)
        if all(abs(y(k:end) - y1) <= band)
            settle_idx = k; break;
        end
    end
    if isnan(settle_idx), settle_time = NaN; else, settle_time = t(settle_idx) - t(k0); end

    M.rise_10_90   = (t90 - t10);
    M.overshoot_frac = overshoot_frac;
    M.t_peak       = t_peak - t(k0);
    M.settle_1pct  = settle_time;
    [settle_time, ~] = settling_time(y, t, k90, y1, A, 1/median(diff(t)));
    M.settle_1pct = settle_time;

end

function k = first_crossing(y, kstart, yth, sgn)
    y = double(y(:));
    k = kstart;
    while k < numel(y) && sgn*(y(k) - yth) < 0
        k = k + 1;
    end
end

function tc = interp_time(t, y, k, yth)
    % Linear interpolation for more precise crossing time
    if k<=1 || k>numel(y), tc = t(min(max(k,1),numel(t))); return; end
    x1=t(k-1); x2=t(k); y1=y(k-1); y2=y(k);
    if y2==y1, tc = x2; else, tc = x1 + (yth - y1)*(x2 - x1)/(y2 - y1); end
end

function [t_settle, k_settle] = settling_time(y, t, k_start, y_final, A, Fs)
    % Noise-aware ±1% band with a minimum dwell time window.
    tol_frac = 0.01;
    W = max(64, round(0.5e-6*Fs));     % need W consecutive in-band samples (≥0.5 µs)
    % estimate noise sigma from the last 10% of the record
    tail = y(round(0.9*numel(y)):end);
    sig  = 1.4826*median(abs(tail - median(tail)));   % MAD→σ
    band = max(tol_frac*abs(A), 4*sig);               % ≥1% OR 4σ (noise-tolerant)

    k = k_start;
    t_settle = NaN; k_settle = NaN;
    while k+W-1 <= numel(y)
        window_ok = all(abs(y(k:k+W-1) - y_final) <= band);
        if window_ok
            t_settle = t(k) - t(k_start); k_settle = k; return;
        end
        k = k + 1;
    end
end

