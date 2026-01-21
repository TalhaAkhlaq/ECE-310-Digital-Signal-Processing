%% Talha Akhlaq
% ECE310 Digital Signal Processing
% Problem Set 1: Spectral Analysis & Sampling
% Questions 2-5
clc; clear; close all;

mag2db_ = @(x) 20*log10(abs(x)+eps);   % magnitude to decibels (prevents log(0))

%% Question 2

normalize_sum1 = @(w) w./sum(w);   % scale window so DC gain = 1

% Rectangular (length 30)
N_rect = 30;
w_rect = normalize_sum1(rectwin(N_rect));

% Chebyshev Window
N_cheb = 95;                  
ripple_dB = 30;                   % sidelobe attenuation (dB)
w_cheb = normalize_sum1(chebwin(N_cheb, ripple_dB));

% Kaiser Window 
N_kais = N_cheb;                  % same as Chebyshev
beta = 3.6;                       % beta controls mainlobe width
w_kais = normalize_sum1(kaiser(N_kais, beta));

[H_rect, w] = freqz(w_rect, 1, 1000);   % rectangular window 
H_cheb = freqz(w_cheb, 1, 1000);        % Chebyshev window 
H_kais = freqz(w_kais, 1, 1000);        % Kaiser window 

figure('Name','Q2: Window Magnitude Responses');
plot(w, mag2db_(H_rect),'LineWidth',1.2); hold on;
plot(w, mag2db_(H_cheb),'LineWidth',1.2);
plot(w, mag2db_(H_kais),'LineWidth',1.2); hold off;
xlim([0, pi]); ylim([-50, 0]); grid on;
xlabel('\omega (rad)'); ylabel('Magnitude (dB)');
legend('Rectangular, N=30','Chebyshev, 30 dB','Kaiser, \beta=3.6','Location','SouthWest');
title('|W(e^{j\omega})| for rectangular, Chebyshev, and Kaiser (DC normalized)');

%% Question 3

% Sum of three continuous-time sinusoids
a   = [1, 0.8, 0.9];                % amplitudes
phi = [pi/2, pi/3, -2*pi/5];        % phases (rad)
f0  = [1e3, 2e3, 3.5e3];            % frequencies (Hz)
fs  = 10e3;                         % sampling rate (Hz)
M   = 1000;                         % number of samples
t   = (0:M-1)/fs;                   % time vector

x = zeros(1,M);
for m = 1:3
    x = x + a(m)*cos(2*pi*f0(m)*t + phi(m));
end

% Hamming window
N = 1024;
w_ham = hamming(M).';
xw = x .* w_ham(1:M);
xwN = zeros(1,N);
xwN(1:min(M,N)) = xw(1:min(M,N));   % adjust to length N
X = fft(xwN, N);

% Bin frequencies 
k    = 0:(N-1);
f_bins = (k*fs)/N;                  % unshifted
Xsh  = fftshift(X);
ksh  = (-N/2):(N/2-1);
f_sh = (ksh*fs)/N;                  % shifted

% Magnitude spectrum (dB)
figure('Name','Q3: Magnitude Spectrum');
plot(f_sh, mag2db_(Xsh),'LineWidth',1.0); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Q3: Magnitude Spectrum (N=1024, Hamming)');

% Nonnegative frequencies only (with Nyquist)
half = 0:(N/2);
f_pos = (half*fs)/N;
X_pos = X(1:(N/2+1));
figure('Name','Q3: Magnitude Spectrum (nonnegative)');
plot(f_pos, mag2db_(X_pos),'LineWidth',1.0); grid on;
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Q3: Magnitude Spectrum, N=1024, Hamming window (nonnegative half)');

% Bin indices (two per tone; positive k)
omega   = 2*pi*f0/fs;               % discrete-time radian freqs
k_ideal = N*f0/fs;                  % bin locations
k_peaks = round([k_ideal, N - k_ideal]);
k_peaks = sort(mod(k_peaks, N));
disp('Q3: Peak DFT bin indices (0..N-1), two per tone (rounded):');
disp(unique(k_peaks));

%% Question 4
Fs = 44100;                            % sampling rate (Hz)
tones = [392, 440, 587.33];            % tone frequencies (Hz)
bpm = 100;                             % tempo (beats per minute)
sec_per_beat = 60/bpm;                 % seconds per quarter note
QN = sec_per_beat;                     % quarter-note duration (s)
samples_per_QN = round(QN*Fs);         % samples per quarter note

% Welch parameters
delta_f = 2;                            % target frequency resolution (Hz)
N0 = round(Fs/delta_f);                 % block length
Nfft = 2^nextpow2(N0);                  
noverlap = floor(N0/2);                 
Mblocks_target = 100;                   % total number of Welch blocks

hop = N0 - noverlap;                      % hop size
L_needed = N0 + (Mblocks_target-1)*hop;   % total samples
num_QN = ceil(L_needed / samples_per_QN); % quarter notes

% Quarter-notes with 2 tones
rng(310);                               
L_total = num_QN * samples_per_QN;
x_audio = zeros(1, L_total);
note_t = (0:samples_per_QN-1)/Fs;
for q = 1:num_QN
    start_idx = (q-1)*samples_per_QN + 1;
    idx = start_idx : start_idx+samples_per_QN-1;
    use = randperm(3,2);
    xn = zeros(1, samples_per_QN);
    for u = use
        xn = xn + cos(2*pi*tones(u)*note_t);   % tones per spec
    end
    x_audio(idx) = xn;
end
% Gaussian noise, SNR=40 dB
signal_power = mean(x_audio.^2);
noise_power = signal_power / (10^(40/10));
x_audio = x_audio + sqrt(noise_power)*randn(size(x_audio));

% 100 Welch blocks                    
x_audio = x_audio(1:L_needed);

% (b) Welch periodogram (Hamming window)
winWelch = hamming(N0);
[Px, fW] = pwelch(x_audio, winWelch, noverlap, Nfft, Fs);  % linear PSD (units/Hz)
Px_dB = 10*log10(Px);                                      % dB/Hz  

% (c) Periodogram Plots
figure('Name','Q4: Welch periodogram');
plot(fW, Px_dB, 'LineWidth', 1.0); grid on;
xlim([0, Fs/2]); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Q4: Welch periodogram (0 to Fs/2)');

figure('Name','Q4: Welch periodogram (300–600 Hz)');
mask = (fW>=300 & fW<=600);
plot(fW(mask), Px_dB(mask), 'LineWidth', 1.0); grid on;
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Q4: Welch periodogram (300–600 Hz)');

% (d) Spectrogram: hop = eighth note
eighth_seconds = QN/2;                  
hop_eighth = round(eighth_seconds*Fs); % samples per hop
hop_spec = min(hop_eighth, N0-1); 
noverlap_spec = max(0, N0 - hop_eighth);
figure('Name','Q4: Spectrogram');
spectrogram(x_audio, winWelch, noverlap_spec, Nfft, Fs);
title('Q4: Spectrogram');
colorbar;

%% Question 5
N = 200;

% Windows
w = ones(1,N);                           % rectangular window
w0 = zeros(1,N);                         % subsampling window
idx_w0 = 4:4:N;                          % select indices
w0(idx_w0) = 1;

% (b) N-point DFTs & Stem Plots
W_N  = fft(w,  N);
W0_N = fft(w0, N);

figure('Name','Q5(b): N-point DFT magnitudes');
subplot(2,1,1);
stem(0:N-1, abs(W_N), 'filled'); grid on;
xlabel('k'); ylabel('|W_N[k]|');
title('Rectangular window, N=200');

subplot(2,1,2);
stem(0:N-1, abs(W0_N), 'filled'); grid on;
xlabel('k'); ylabel('|W0_N[k]|');
title('w0[n], N=200');

% (c) NO-point DFTs
N0_hr = 2^nextpow2(16*N);               % smallest power
W_hr  = fft(w,  N0_hr);
W0_hr = fft(w0, N0_hr);
W_hr  = W_hr  / W_hr(1);                % normalize DC to 1
W0_hr = W0_hr / W0_hr(1);

% Frequency grid 
k0 = 0:N0_hr-1;
omega0 = 2*pi*k0/N0_hr;                 % unshifted
omega_sh = omega0 - 2*pi*(omega0>pi);   % shifted
W_hr_sh  = fftshift(W_hr);
W0_hr_sh = fftshift(W0_hr);

figure('Name','Q5(c): Magnitude spectra');
plot(omega_sh, abs(W_hr_sh),'LineWidth',1.1); hold on;
plot(omega_sh, abs(W0_hr_sh),'LineWidth',1.1); hold off; grid on;
xlabel('\omega (rad)'); ylabel('Magnitude');
legend('|W(\omega)|','|W_0(\omega)|','Location','Best');
title('Q5(c): |W(\omega)| and |W_0(\omega)|, −\pi ≤ \omega ≤ \pi');


% (d) N0-point DFTs
A = 1; theta = 0; omega0_sig = 0.4;
nN = 0:N-1;
x = A*cos(omega0_sig*nN + theta);

xw_rect = x .* w;                      % rectangular window
xw_w0   = x .* w0;                     % subsampling window
X_rect = fft(xw_rect, N0_hr);          % N0-point DFT (rectangular)
X_w0   = fft(xw_w0,   N0_hr);          % N0-point DFT (w0)

% Plot magnitudes
mask_half = omega0 <= pi + 1e-12;      
figure('Name','Q5(d): Magnitude spectra');
plot(omega0(mask_half), abs(X_rect(mask_half)),'LineWidth',1.0); hold on;
plot(omega0(mask_half), abs(X_w0(mask_half)),'LineWidth',1.0); hold off; grid on;
xlabel('\omega (rad)'); ylabel('Magnitude');
legend('Rectangular window','w0[n] window','Location','Best');
title('Q5(d): |X(\omega)| for 0 ≤ \omega ≤ \pi');

% (e) Comparison
disp('Q5(e): w0[n] introduces imaging distortion from subsampling. Spectral resolution is not qualitatively different compared to the rectangular window')
