%% Talha Akhlaq
% ECE310 Digital Signal Processing
% Problem Set 2: Filters
clear; close all; clc;

%% Specifications
Fs = 40e6; % Hz
f_pass = [10e6 11e6]; % Hz
f_stop = [9e6 12e6]; % Hz
Rp = 1.5; % dB passband variation
Rs = 30; % dB stopband attenuation
NYQ = Fs/2; % Hz
f_axis_MHz = linspace(0, 20, 4001); % MHz axis for plotting & sampling responses
f_axis_Hz  = 1e6 * f_axis_MHz; % Hz
w_axis_rad = 2*pi*f_axis_Hz; % rad/s (analog)

Wp_dig = f_pass/NYQ; % normalized to Nyquist for digital IIR
Ws_dig = f_stop/NYQ;
Wp_an  = 2*pi*f_pass; % rad/s for analog IIR
Ws_an  = 2*pi*f_stop; % rad/s

%% dB and deg 
mag2db_ = @(x) 20*log10(abs(x));  
rad2deg_= @(x) (180/pi)*x;

%% Question 1
iir_kinds = {'butter','cheby1','cheby2','ellip'};
Q1_results = struct(); % collect metrics

% Pre-allocate
all_mag_dB = []; all_phase_deg = [];

% Determine global limits
resp_cache = struct();

for kd = 1:2 % 1=analog, 2=digital
    dom = ternary(kd==1,'Analog','Digital');
    for kt = 1:numel(iir_kinds)
        kind = iir_kinds{kt};

        if kd==1 % Analog
            switch kind
                case 'butter'
                    [n,Wn] = buttord(Wp_an, Ws_an, Rp, Rs, 's');
                    [z,p,k] = butter(n, Wn, 'bandpass', 's');
                case 'cheby1'
                    [n,Wn] = cheb1ord(Wp_an, Ws_an, Rp, Rs, 's');
                    [z,p,k] = cheby1(n, Rp, Wn, 'bandpass', 's');
                case 'cheby2'
                    [n,Wn] = cheb2ord(Wp_an, Ws_an, Rp, Rs, 's');
                    [z,p,k] = cheby2(n, Rs, Wn, 'bandpass', 's');
                case 'ellip'
                    [n,Wn] = ellipord(Wp_an, Ws_an, Rp, Rs, 's');
                    [z,p,k] = ellip(n, Rp, Rs, Wn, 'bandpass', 's');
            end
            % Transfer function
            [b,a] = zp2tf(z,p,k);
            H = freqs(b,a,w_axis_rad);
        else % Digital
            switch kind
                case 'butter'
                    [n,Wn] = buttord(Wp_dig, Ws_dig, Rp, Rs);
                    [z,p,k] = butter(n+1, Wn, 'bandpass');
                case 'cheby1'
                    [n,Wn] = cheb1ord(Wp_dig, Ws_dig, Rp, Rs);
                    Rp_eff = 1.49;
                    [z,p,k] = cheby1(n, Rp_eff, Wn, 'bandpass');
                case 'cheby2'
                    [n,Wn] = cheb2ord(Wp_dig, Ws_dig, Rp, 30.1);
                    [z,p,k] = cheby2(n+4, 30.1, Wn, 'bandpass'); 
                case 'ellip'
                    [n,Wn] = ellipord(Wp_dig, Ws_dig, Rp, Rs);
                    Rp_eff = 1.49;
                    [z,p,k] = ellip(n, Rp_eff, Rs, Wn, 'bandpass'); 
            end
            [b,a] = zp2tf(z,p,k);
            Nw = 8192;
            [Hf, ff] = freqz(b,a,Nw,Fs); % ff in Hz
            H = interp1(ff, Hf, f_axis_Hz, 'linear', 'extrap'); % map to 0..20 MHz
        end

        mag_dB  = mag2db_(H);
        ph_deg  = rad2deg_(unwrap(angle(H)));

        all_mag_dB  = [all_mag_dB; mag_dB(:).'];
        all_phase_deg = [all_phase_deg; ph_deg(:).'];

        key = sprintf('%s_%s', lower(dom), lower(kind));
        resp_cache.(key).z = z;
        resp_cache.(key).p = p;
        resp_cache.(key).k = k;
        resp_cache.(key).b = b;
        resp_cache.(key).a = a;
        resp_cache.(key).n = n;
        resp_cache.(key).mag_dB = mag_dB;
        resp_cache.(key).ph_deg = ph_deg;
        resp_cache.(key).dom = dom;
        resp_cache.(key).kind = kind;
    end
end

% Determine y-limits
mag_min = min(all_mag_dB, [], 'all'); mag_max = max(all_mag_dB, [], 'all');
% Cap magnitude upper limit to just below 0 dB for better comparison
mag_ylim = [min(-100, floor(mag_min/5)*5), 5];
% Set phase y-limits from global min/max with small padding
phi_min = min(all_phase_deg, [], 'all'); phi_max = max(all_phase_deg, [], 'all');
pad = 0.05*(phi_max - phi_min + eps);
phi_ylim = [phi_min - pad, phi_max + pad];

% Plotting
for kd = 1:2
    dom = ternary(kd==1,'Analog','Digital');
    for kt = 1:numel(iir_kinds)
        kind = iir_kinds{kt};
        key = sprintf('%s_%s', lower(dom), lower(kind));
        z = resp_cache.(key).z; p = resp_cache.(key).p; k = resp_cache.(key).k;
        b = resp_cache.(key).b; a = resp_cache.(key).a; n = resp_cache.(key).n;
        mag_dB = resp_cache.(key).mag_dB; ph_deg = resp_cache.(key).ph_deg;

        % Question 1 (a)
        order_den = length(a)-1;

        % Question 1 (b)
        figure('Name', sprintf('Q1 %s %s — PZ', dom, kind));
        if kd==1
            % Analog: s-plane plot using z,p; axes in rad/s
            plot(real(p), imag(p), 'x', 'MarkerSize', 8); hold on;
            if ~isempty(z), plot(real(z), imag(z), 'o', 'MarkerSize', 6); end
            grid on; axis equal;
            xlabel('Real(s) [rad/s]'); ylabel('Imag(s) [rad/s]');
            title(sprintf('Analog %s — Poles (x) and Zeros (o)', upper(kind)));
        else
            % Digital: z-plane with coefficients 
            zplane(b, a); grid on;
            title(sprintf('Digital %s — Pole-Zero Plot (Unit Circle)', upper(kind)));
        end

        % Question 1 (c)
        figure('Name', sprintf('Q1 %s %s — Resp', dom, kind));
        subplot(2,1,1);
        plot(f_axis_MHz, mag_dB, 'LineWidth', 1.0); grid on;
        xlabel('Frequency [MHz]'); ylabel('Magnitude [dB]');
        title(sprintf('%s %s: |H| (dB)', dom, upper(kind)));
        ylim(mag_ylim);
        xlim([0 20]);
        subplot(2,1,2);
        plot(f_axis_MHz, ph_deg, 'LineWidth', 1.0); grid on;
        xlabel('Frequency [MHz]'); ylabel('Unwrapped phase [deg]');
        title(sprintf('%s %s: H (deg)', dom, upper(kind)));
        ylim(phi_ylim);
        xlim([0 20]);

        % Question 1 (d)
        idx_pb = f_axis_Hz >= f_pass(1) & f_axis_Hz <= f_pass(2);
        idx_sbL= f_axis_Hz <= f_stop(1); % Unused
        idx_sbU= f_axis_Hz >= f_stop(2); % Unused
        [pk_pb_dB,~] = max(mag_dB(idx_pb));
        mag_dB_at = @(fMHz) interp1(f_axis_MHz, mag_dB, fMHz, 'linear');
        pb10 = mag_dB_at(10); pb11 = mag_dB_at(11);
        sb9  = mag_dB_at(9);  sb12 = mag_dB_at(12);

        meets_pb = (max(abs([pb10 pb11])) <= Rp+1e-3);
        meets_sb = (sb9  <= -Rs+1e-6) && (sb12 <= -Rs+1e-6);
        overdesign_sb = (sb9 < -Rs) && (sb12 < -Rs);

        Q1_results.(key) = struct('domain',dom,'kind',kind,'n',n,'order_den',order_den, ...
            'pk_pb_dB',pk_pb_dB,'pb10_dB',pb10,'pb11_dB',pb11, ...
            'sb9_dB',sb9,'sb12_dB',sb12,'meets_pb',meets_pb, ...
            'meets_sb',meets_sb,'overdesign_sb',overdesign_sb);

        % Results: peak passband gain ~ 0 dB
        % Edges (10, 11, 9, 12 MHz): all 8 satisfy Rp=1.5 dB and Rs=30 dB
        % Overdesign: Butterworth and Chebyshev I stopbands by ~4-8 dB
        %             Analog Chebyshev II large on one stopband (zeros)
        %             Elliptic and digital Chebyshev II sit near 30 dB
        % Cause: Butterworth monotone; Chebyshev I passband equiripple; Chebyshev II stopband zeros; Elliptic equiripple both bands
        % Tweaks: digital Butterworth n+1; Chebyshev I Rp=1.49 dB; Chebyshev II Rs=30.1 dB and n+4; Elliptic Rp=1.49 dB
    end
end

% Output
kinds   = {'butter','cheby1','cheby2','ellip'};
domains = {'Analog','Digital'};

for kd = 1:2
    dom = domains{kd};
    fprintf('\nQuestion 1 %s:\n', dom);
    fprintf('  %-8s  %5s  %6s  %9s  %9s  %9s  %9s  %7s  %7s  %9s\n', ...
        'Type','Npro','Nreal','pb@10M','pb@11M','sb@9M','sb@12M','PB_ok','SB_ok','pkPB');

    for kt = 1:numel(kinds)
        kind = kinds{kt};
        key  = sprintf('%s_%s', lower(dom), kind);
        R    = Q1_results.(key);

        fprintf('  %-8s  %5d  %6d  %9.4f  %9.4f  %9.4f  %9.4f  %7d  %7d  %9.4f\n', ...
            upper(kind), R.n, R.order_den, R.pb10_dB, R.pb11_dB, ...
            R.sb9_dB, R.sb12_dB, R.meets_pb, R.meets_sb, R.pk_pb_dB);
    end
end

fprintf('\n\n')

%% Question 2
% Convert dB specs to linear FIR
Ap = 10^(Rp/20);                 
dp = (Ap - 1)/(Ap + 1);          
ds = 10^(-Rs/20);                

% Bands: [0 fstop1] [fpass1 fpass2] [fstop2 Nyq]
F_edges = [0, f_stop(1), f_pass(1), f_pass(2), f_stop(2), NYQ]; % Unused
A_bands = [0, 1, 0];
Dev     = [ds, dp, ds];

% Question 2 (a)
[nK, WnK, betaK, ftypeK] = kaiserord([f_stop(1) f_pass(1) f_pass(2) f_stop(2)], A_bands, Dev, Fs);
nK = max(nK, 2*ceil(nK/2)); % even order for symmetry
for iter = 1:100
    hK = fir1(nK, WnK, ftypeK, kaiser(nK+1, betaK), 'scale');
    [Hk_d, fk_d] = freqz(hK, 1, 65536, Fs);
    Mk_dB = 20*log10(abs(Hk_d));
    in_pb  = fk_d >= f_pass(1) & fk_d <= f_pass(2);
    in_sbL = fk_d <= f_stop(1); 
    in_sbU = fk_d >= f_stop(2);
    pb_spread = max(Mk_dB(in_pb)) - min(Mk_dB(in_pb));
    sbL_peak  = max(Mk_dB(in_sbL));
    sbU_peak  = max(Mk_dB(in_sbU));
    if pb_spread <= Rp && sbL_peak <= -Rs && sbU_peak <= -Rs
        break
    end
    nK = nK + 2;
end

% Pole-zero via zplane; stem of coefficients; frequency response
figure('Name','Q2 Kaiser — stem(h)'); stem(0:length(hK)-1, hK, 'filled'); grid on;
xlabel('n'); ylabel('h[n]'); title(sprintf('Kaiser FIR, length %d', length(hK)));

figure('Name','Q2 Kaiser — PZ'); zplane(hK, 1); grid on;
title('Kaiser FIR — Pole-Zero Plot');

[Hk_full, fk_full] = freqz(hK, 1, 8192, Fs);
Hk = interp1(fk_full, Hk_full, f_axis_Hz, 'linear','extrap');
magK_dB = mag2db_(Hk); phK_deg = rad2deg_(unwrap(angle(Hk)));

figure('Name','Q2 Kaiser — Resp');
subplot(2,1,1); plot(f_axis_MHz, magK_dB, 'LineWidth',1.0); grid on;
xlabel('Frequency [MHz]'); ylabel('Magnitude [dB]');
title('Kaiser FIR: |H| (dB)'); ylim(mag_ylim); xlim([0 20]);
subplot(2,1,2); plot(f_axis_MHz, phK_deg, 'LineWidth',1.0); grid on;
xlabel('Frequency [MHz]'); ylabel('Unwrapped phase [deg]');
title('Kaiser FIR: H (deg)'); ylim(phi_ylim); xlim([0 20]);

% Equiripple FIR
[nE, foE, aoE, wE] = firpmord([f_stop(1) f_pass(1) f_pass(2) f_stop(2)], A_bands, Dev, Fs);
if numel(wE) ~= 3, wE = 1./[ds dp ds]; end
nE = max(nE, 2*ceil(nE/2));
for iter = 1:100
    bE = firpm(nE, foE, aoE, wE);
    [He_d, fe_d] = freqz(bE, 1, 65536, Fs);
    Me_dB = 20*log10(abs(He_d));
    in_pb  = fe_d >= f_pass(1) & fe_d <= f_pass(2);
    in_sbL = fe_d <= f_stop(1);
    in_sbU = fe_d >= f_stop(2);
    pb_spread = max(Me_dB(in_pb)) - min(Me_dB(in_pb));
    sbL_peak  = max(Me_dB(in_sbL));
    sbU_peak  = max(Me_dB(in_sbU));
    if pb_spread <= Rp && sbL_peak <= -Rs && sbU_peak <= -Rs
        break
    end
    nE = nE + 2;
end

figure('Name','Q2 Equiripple — stem(h)'); stem(0:length(bE)-1, bE, 'filled'); grid on;
xlabel('n'); ylabel('h[n]'); title(sprintf('Equiripple FIR, length %d', length(bE)));

figure('Name','Q2 Equiripple — PZ'); zplane(bE, 1); grid on;
title('Equiripple FIR — Pole-Zero Plot');

[He_full, fe_full] = freqz(bE, 1, 8192, Fs); 
He = interp1(fe_full, He_full, f_axis_Hz, 'linear','extrap');
magE_dB = mag2db_(He); phE_deg = rad2deg_(unwrap(angle(He)));

figure('Name','Q2 Equiripple — Resp');
subplot(2,1,1); plot(f_axis_MHz, magE_dB, 'LineWidth',1.0); grid on;
xlabel('Frequency [MHz]'); ylabel('Magnitude [dB]');
title('Equiripple FIR: |H| (dB)'); ylim(mag_ylim); xlim([0 20]);
subplot(2,1,2); plot(f_axis_MHz, phE_deg, 'LineWidth',1.0); grid on;
xlabel('Frequency [MHz]'); ylabel('Unwrapped phase [deg]');
title('Equiripple FIR: H (deg)'); ylim(phi_ylim); xlim([0 20]);

% Question 2(b)
wE_full = wE(:);
if numel(wE_full) ~= 3, wE_full = (1./Dev(:)); end
weights_vs_tolerances = table( ...
    {'Stopband L';'Passband';'Stopband U'}, ...
    Dev(:), ...
    wE_full, ...
    'VariableNames', {'Band','Deviation_linear','Weight_used'});
disp('Question 2(b):'); 
disp(weights_vs_tolerances);

% Question 2(c)
check_metrics = @(mag_dB) struct( ...
    'pb_max_dB', max(mag_dB(f_axis_Hz >= f_pass(1) & f_axis_Hz <= f_pass(2))), ...
    'pb_min_dB', min(mag_dB(f_axis_Hz >= f_pass(1) & f_axis_Hz <= f_pass(2))), ...
    'pb_spread_dB', max(mag_dB(f_axis_Hz >= f_pass(1) & f_axis_Hz <= f_pass(2))) - ...
                   min(mag_dB(f_axis_Hz >= f_pass(1) & f_axis_Hz <= f_pass(2))), ...
    'sbL_peak_dB', max(mag_dB(f_axis_Hz <= f_stop(1))), ...
    'sbU_peak_dB', max(mag_dB(f_axis_Hz >= f_stop(2))) );

metK = check_metrics(magK_dB);
metE = check_metrics(magE_dB);
disp('Question 2(c) Kaiser:'); disp(metK);
disp('Question 2(c) Equiripple:'); disp(metE);

% Comments:
% Passband: Kaiser 0.223 dB (1.277 dB under 1.5); Equiripple 1.448 dB (0.052 dB under).
% Stopbands: Both below −30 dB.
% Both meet specs; Kaiser is overdesigned

%% ternary
function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end
