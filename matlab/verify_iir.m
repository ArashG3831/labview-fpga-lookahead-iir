% ================= Timestamp-only check (no post-fit/no xcorr) =================
% Goal: Prove FPGA filter output matches MATLAB reference with ONLY:
%   - timestamp alignment (same Fs),
%   - known fixed pipeline delay of -3 samples.
% No cross-correlation, no sign flip, no gain/offset fitting.

clear; clc; close all;

% ----------- Inputs (two columns: time, value) -----------
raw_file  = 'raw2.txt';   % pre-filter ADC stream (t, x)
fpga_file = 'fi2.txt';    % FPGA-filtered stream (t, y_fpga)

% Your look-ahead biquad (k=1)
b = [0.04613180209331293, 0.15257101841705373,...
     0.18943481534878998, 0.10568378381967050, 0.02268818479462131];
a = [1, 0, -0.7253696722084265, 0, 0.24187927668187462];

% Known deterministic pipeline delay on FPGA output:
PIPE_DELAY = -3;  % samples (negative means FPGA lags ideal by 3 samples)

% Analysis settings (just to avoid edge artifacts in plots/metrics)
EDGE_TRIM  =  max(32, 1024);     % drop this many samples at each end
NOISE_BAND = [1e6, 60e6];        % Hz

% ---------------------- Load & stitch (timestamp only) -------------------
[t,  x_raw,  Fs]  = load_trace_stitched(raw_file);
[tf, y_hw,   FsF] = load_trace_stitched(fpga_file);
assert(abs(Fs-FsF) <= 1e-9*Fs, 'Fs mismatch: raw=%.9g, fpga=%.9g', Fs, FsF);

N  = numel(t);
fprintf('Fs = %.6f MHz, N=%d\n', Fs/1e6, N);

% Align streams by integer offset from their start timestamps (no xcorr)
start_offset_samples = round((tf(1) - t(1)) * Fs);
y_on_raw_grid        = circshift(y_hw, -start_offset_samples);

% Apply ONLY the known fixed pipeline delay
y_fpga_aligned = circshift(y_on_raw_grid, PIPE_DELAY);

% ---------------------- MATLAB reference (zero ICs) ----------------------
% Same initial conditions: start MATLAB IIR from zero state
y_ideal = filter(b, a, double(x_raw));

% ---------------------- Window for metrics (steady-state) ----------------
edge = min(EDGE_TRIM, floor(N/8));
win  = (1+edge):(N-edge);

% ---------------------- Time plot (no scaling) ---------------------------
tt = (t - t(1))*1e6;
figure('Color','w','Name','Timestamp-only alignment (delay = -3)');
plot(tt(win), y_ideal(win),        'k','LineWidth',1.1); hold on;
plot(tt(win), y_fpga_aligned(win), 'r--','LineWidth',1.0);
grid on; xlabel('Time (\mus)'); ylabel('Amplitude (a.u.)');
title(sprintf('Filters only — FPGA (unaltered) vs MATLAB  [delay=%+d]', PIPE_DELAY));
legend('MATLAB ideal','FPGA (timestamp + fixed delay)','Location','best');

% ---------------------- Error & metrics (no gain/offset) -----------------
e        = y_fpga_aligned(win) - y_ideal(win);
rel_rms  = rms(e) / rms(y_ideal(win));
nmse_db  = 20*log10(rel_rms);
fprintf('Relative RMS error (no post-fit) = %.4g  (NMSE = %.2f dB)\n', rel_rms, nmse_db);

% ---------------------- PSDs (one-sided Welch) ---------------------------
[f, Praw] = welch_psd(x_raw(win),         Fs);
[~, Pid ] = welch_psd(y_ideal(win),       Fs);
[~, Pfp ] = welch_psd(y_fpga_aligned(win),Fs);

figure('Color','w','Name','Welch PSD (no post-fit)');
plot(f/1e6, 10*log10(Praw+realmin),'b','LineWidth',1.0); hold on;
plot(f/1e6, 10*log10(Pfp +realmin),'r','LineWidth',1.15);
plot(f/1e6, 10*log10(Pid +realmin),'k','LineWidth',1.15);
grid on; xlim([f(2) Fs/2]/1e6);
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)');
title('Welch PSD — Raw (blue) vs FPGA (red) vs MATLAB ideal (black)');
legend('Raw','FPGA (timestamp+delay only)','MATLAB ideal','Location','best');

% Band-integrated noise (identical method for both)
mask = (f >= NOISE_BAND(1) & f <= NOISE_BAND(2));
Vr = sqrt(trapz(f(mask), Praw(mask)));
Vf = sqrt(trapz(f(mask), Pfp (mask)));
Vm = sqrt(trapz(f(mask), Pid (mask)));
fprintf('Band %.1f–%.1f MHz: Raw=%.3g | FPGA=%.3g | MATLAB=%.3g  (no post-fit)\n', ...
    NOISE_BAND(1)/1e6, NOISE_BAND(2)/1e6, Vr, Vf, Vm);

% Optional: Error PSD (diagnostic)
[fe, Pe] = welch_psd(e, Fs);
figure('Color','w','Name','Error PSD (no post-fit)');
plot(fe/1e6, 10*log10(Pe+realmin),'LineWidth',1.1); grid on;
xlabel('Frequency (MHz)'); ylabel('PSD (dB/Hz)'); title('Error PSD');

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

function [f, P] = welch_psd(x, Fs)
    x = x(:) - mean(x);
    N  = numel(x);
    L  = max(256, 2^floor(log2(max(N/4,256))));
    ov = floor(0.5*L);
    nF = max(2048, 2^nextpow2(L));
    [P, f] = pwelch(x, hann(L,'periodic'), ov, nF, Fs, 'onesided');
end
