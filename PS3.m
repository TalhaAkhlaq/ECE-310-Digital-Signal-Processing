%% Talha Akhlaq
% ECE310 Digital Signal Processing
% Problem Set 3
clear; close all; clc;
set(groot,'defaultAxesXGrid','on');
set(groot,'defaultAxesYGrid','on');
set(groot,'defaultAxesXMinorGrid','on');
set(groot,'defaultAxesYMinorGrid','on');

%% Problem 4(a)

% Parameters
Rp = 1.5;
Rs = 30;
Wp = [0.3 0.6];
n = 4;

% Discrete Time
[z,p,k] = ellip(n, Rp, Rs, Wp, 'bandpass');
[b,a] = zp2tf(z,p,k);

N = max(numel(a), numel(b)) - 1;
fprintf('4(a) Bandpass Order: %d\n', N);

% Magnitude Response
nfft = 4096;
[H,w] = freqz(b, a, nfft);
figure;
plot(w/pi, 20*log10(abs(H)), 'LineWidth', 1.2, 'Color', [0.10 0.75 0.65]);
grid on;
xlabel('\omega/\pi');
ylabel('Magnitude (dB)');
title('4(a): Elliptic Bandpass, N=8');
ylim([-110 5]);
xlim([0 1]);
hold on;
yl = ylim;
plot([Wp;Wp],[yl;yl],'k--');
hold off;

%% Problem 4(b)

[sos_up,   g_up]   = zp2sos(z, p, k, 'up',   'inf');
[sos_down, g_down] = zp2sos(z, p, k, 'down', 'inf');

[bu, au] = sos2tf(sos_up,   g_up);
[bd, ad] = sos2tf(sos_down, g_down);
fprintf('\n');
fprintf('4(b) Up:   max|change in b| = %g,  max|change in a| = %g\n', max(abs(bu-b)), max(abs(au-a)));
fprintf('4(b) Down: max|change in b| = %g,  max|change in a| = %g\n', max(abs(bd-b)), max(abs(ad-a)));

%% Problem 4(c)

sosU = sos_up;
sosD = sos_down;
L = size(sosU,1);

stage_mags = @(row) deal(sort(abs(roots(row(1:3))).','ascend'), sort(abs(roots(row(4:6))).','ascend'));

% Per-Stage Magnitudes
mzU = zeros(L,2);
mpU = zeros(L,2);
mzD = zeros(L,2);
mpD = zeros(L,2);
for i = 1:L
    [mzU(i,:), mpU(i,:)] = stage_mags(sosU(i,:));
    [mzD(i,:), mpD(i,:)] = stage_mags(sosD(i,:));
end

% Pole Magnitudes
fprintf('\n');
fprintf('4(c) Pole Magnitudes (Up):\n');
for i = 1:L
    fprintf('%.6f %.6f\n', mpU(i,1), mpU(i,2));
end

fprintf('\n4(c) Pole Magnitudes (Down):\n');
for i = 1:L
    fprintf('%.6f %.6f\n', mpD(i,1), mpD(i,2));
end

% Zero Magnitudes
fprintf('\n4(c) Zero Magnitudes (Up):\n');
for i = 1:L
    fprintf('%.6f %.6f\n', mzU(i,1), mzU(i,2));
end

fprintf('\n4(c) Zero Magnitudes (Down):\n');
for i = 1:L
    fprintf('%.6f %.6f\n', mzD(i,1), mzD(i,2));
end

% Monotonicity Checks
tol = 1e-3;
mup = all(diff(mean(mpU,2)) >= -tol);
mdp = all(diff(mean(mpD,2)) <= tol);
muz = all(diff(mean(mzU,2)) >= -tol);
mdz = all(diff(mean(mzD,2)) <= tol);
fprintf('\n4(c) Monotonicity Checks:\n');
fprintf('Up Poles Increasing   = %d\n',   mup);
fprintf('Down Poles Decreasing = %d\n', mdp);
fprintf('Up Zeros Increasing   = %d\n',   muz);
fprintf('Down Zeros Decreasing = %d\n', mdz);

%% Problem 4(d)

nfft = 8192;
m_up   = plotFkCumulative(sos_up, g_up, nfft, '4(d): Cumulative F_k (Up Order)', [-130 10]);
fprintf('\n4(d) max |F_k| (Up):\n');
fprintf('%.4f\n', m_up);
fprintf('\n');

m_down = plotFkCumulative(sos_down, g_down, nfft, '4(d): Cumulative F_k (Down Order)', [-120 10]);
fprintf('4(d) max |F_k| (Down):\n');
fprintf('%.4f\n', m_down);
fprintf('\n');

for f = 1:3
    figure(f);
    ax = gca;
    yl = ax.YLim;
    ax.YTick = ceil(yl(1)/10)*10 : 10 : floor(yl(2)/10)*10;
end

%% Problem 4(e)

Bup = sos_up(:,1:3);
Aup = sos_up(:,4:6);
Bdn = sos_down(:,1:3);
Adn = sos_down(:,4:6);

L = size(Aup,1);
fprintf('4(e) Up vs Down SOS (Reverse Order):\n\n');

for k = 1:L
    j = L-k+1;

    % Denominator
    den_err = norm(Aup(k,:) - Adn(j,:), 2) / max(norm(Aup(k,:),2), eps);

    % Least-Squares Scaling Factor
    d = Bdn(j,:).';
    u = Bup(k,:).';
    c = (d'*u) / (d'*d);
    num_err = norm(Bup(k,:) - real(c)*Bdn(j,:), 2) / max(norm(Bup(k,:),2), eps);

    fprintf('Up %d vs Down %d\n', k, j);
    fprintf('  Denominator Error: %.3e\n', den_err);
    fprintf('  Scale c: %.5g\n', real(c));
    fprintf('  Numerator Error: %.3e\n\n', num_err);
end

%% Problem 5(a)

bH = [0.1336 0.0563 0.0563 0.1336];
aH = [1.0000 -1.5055 1.2630 -0.3778];
b1 = [-0.4954 1.0000];
a1 = [1.0000 -0.4954];
b2 = [ 0.7626 -1.0101 1.0000];
a2 = [1.0000 -1.0101 0.7626];
bHa = 0.5*(conv(b1,a2) + conv(b2,a1));
aHa = conv(a1,a2);

% Frequency Response Check
nfft = 4096;
[H1,w] = freqz(bH, aH, nfft);
[Ha,w] = freqz(bHa,aHa,nfft);
max_resp_err = max(abs(H1 - Ha));
fprintf('5(a) Max |H - H_a| = %.6e\n', max_resp_err);

% Pole Comparison
zH = roots(bH);
pH = roots(aH);
zHa = roots(bHa);
pHa = roots(aHa);

DhZ = abs(zH - zHa.');
DhP = abs(pH - pHa.');
max_zero_err = max([max(min(DhZ,[],2)) max(min(DhZ,[],1))]);
max_pole_err = max([max(min(DhP,[],2)) max(min(DhP,[],1))]);

fprintf('5(a) Max Zero Error = %.6e\n', max_zero_err);
fprintf('5(a) Max Pole Error = %.6e\n', max_pole_err);

%% Problem 5(b)

% 4-bit Two's Complement
w = 4;

% H(z) Coefficients
bH = [0.1336 0.0563 0.0563 0.1336];
aH = [1.0000 -1.5055 1.2630 -0.3778];
fH = choose_frac_len([bH aH], w);
bHQ = quantize_twos(bH, w, fH);
aHQ = quantize_twos(aH, w, fH);
fprintf('\n5(b) H(z) Fractional Bits = %d\n', fH);

% Allpass Pair Coefficients
b1 = [-0.4954 1.0000];
a1 = [1.0000 -0.4954];
b2 = [ 0.7626 -1.0101 1.0000];
a2 = [1.0000 -1.0101 0.7626];

f1 = choose_frac_len([b1 a1], w);
f2 = choose_frac_len([b2 a2], w);
b1q = quantize_twos(b1, w, f1);
a1q = quantize_twos(a1, w, f1);
b2q = quantize_twos(b2, w, f2);
a2q = quantize_twos(a2, w, f2);
fprintf('5(b) Allpass1 Fractional Bits = %d\n', f1);
fprintf('5(b) Allpass2 Fractional Bits = %d\n\n', f2);

% Quantized Composites
bHaQ = 0.5*(conv(b1q,a2q) + conv(b2q,a1q));
aHaQ = conv(a1q,a2q);

%% Problem 5(c)

% Poles
zH = roots(bH);
pH = roots(aH);
zHQ = roots(bHQ);
pHQ = roots(aHQ);
zHaQ = roots(bHaQ);
pHaQ = roots(aHaQ);

% X Plane Overlay
figure;
hold on;
th = linspace(0,2*pi,512);
plot(cos(th), sin(th), 'w--', 'LineWidth', 1.0);
plot(real(zH), imag(zH), 'bo', 'LineWidth', 1.2);
plot(real(pH), imag(pH), 'bx', 'LineWidth', 1.2);
plot(real(zHQ), imag(zHQ), 'ro', 'LineWidth', 1.2);
plot(real(pHQ), imag(pHQ), 'rx', 'LineWidth', 1.2);
plot(real(zHaQ), imag(zHaQ), 'go', 'LineWidth', 1.2);
plot(real(pHaQ), imag(pHaQ), 'gx', 'LineWidth', 1.2);
grid on;
axis equal;
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
xlabel('Real Part');
ylabel('Imaginary Part');
title('5(c): Pole-Zero Map');
legend('Unit circle','H zeros','H poles','H_Q zeros','H_Q poles','H_{aQ} zeros','H_{aQ} poles','Location','best');

% Maximum Absolute Errors
hausdorff = @(u,v) max([max(min(abs(u - v.'),[],2)) max(min(abs(u - v.'),[],1))]);

max_pole_err_H_HQ = hausdorff(pH, pHQ);
max_pole_err_H_HaQ = hausdorff(pH, pHaQ);
max_zero_err_H_HQ = hausdorff(zH, zHQ);
max_zero_err_H_HaQ = hausdorff(zH, zHaQ);

fprintf('5(c) Max Pole Error (H vs H_Q)  = %.3e\n', max_pole_err_H_HQ);
fprintf('5(c) Max Pole Error (H vs H_{aQ}) = %.3e\n', max_pole_err_H_HaQ);
fprintf('5(c) Max Zero Error (H vs H_Q)  = %.3e\n', max_zero_err_H_HQ);
fprintf('5(c) Max Zero Error (H vs H_{aQ}) = %.3e\n\n', max_zero_err_H_HaQ);

%% Problem 5(d)

% Gain
evaltf = @(b,a,z) polyval(b,z)./polyval(a,z);
todB = @(x) 20*log10(max(abs(x),1e-12));

% Original
H0 = evaltf(bH, aH, 1);
Hpi = evaltf(bH, aH, -1);
H0dB = todB(H0);
HpidB = todB(Hpi);

% Quantized H_Q
HQ0 = evaltf(bHQ, aHQ, 1);
HQpi = evaltf(bHQ, aHQ, -1);
HQ0dB = todB(HQ0);
HQpidB = todB(HQpi);

% Quantized H_aQ
HaQ0 = evaltf(bHaQ, aHaQ, 1);
HaQpi = evaltf(bHaQ, aHaQ, -1);
HaQ0dB = todB(HaQ0);
HaQpidB = todB(HaQpi);

% Errors
err0_HQ = HQ0dB - H0dB;
errpi_HQ = HQpidB - HpidB;
err0_HaQ = HaQ0dB - H0dB;
errpi_HaQ = HaQpidB - HpidB;

fprintf('5(d) Original Gains:\n');
fprintf('  |H(0)|  = %.6f (%.3f dB)\n', abs(H0), H0dB);
fprintf('  |H(pi)| = %.6e (%.3f dB)\n\n', abs(Hpi), HpidB);

fprintf('5(d) H_Q Gains:\n');
fprintf('  |H_Q(0)|  = %.6f (%.3f dB)\n', abs(HQ0), HQ0dB);
fprintf('  |H_Q(pi)| = %.6e (%.3f dB)\n\n', abs(HQpi), HQpidB);

fprintf('5(d) H_{aQ} Gains:\n');
fprintf('  |H_{aQ}(0)|  = %.6f (%.3f dB)\n', abs(HaQ0), HaQ0dB);
fprintf('  |H_{aQ}(pi)| = %.6e (%.3f dB)\n\n', abs(HaQpi), HaQpidB);

fprintf('5(d) dB Error (H_Q vs H):\n');
fprintf('  omega=0:  %.3e dB\n', err0_HQ);
fprintf('  omega=pi: %.3e dB\n\n', errpi_HQ);

fprintf('5(d) dB Error (H_{aQ} vs H):\n');
fprintf('  omega=0:  %.3e dB\n', err0_HaQ);
fprintf('  omega=pi: %.3e dB\n\n', errpi_HaQ);

%% Problem 5(e)

nfft = 10000;
[H,w] = freqz(bH, aH, nfft);
[HQ,~] = freqz(bHQ, aHQ, nfft);
[HaQ,~] = freqz(bHaQ, aHaQ, nfft);

max_err_H_HQ = max(abs(H - HQ));
max_err_H_HaQ = max(abs(H - HaQ));

fprintf('5(e) max|H - H_Q|    = %.3e\n', max_err_H_HQ);
fprintf('5(e) max|H - H_{aQ}| = %.3e\n\n', max_err_H_HaQ);

%% Problem 5(f)

nfft = 10000;
[H,w] = freqz(bH, aH, nfft);
[HQ,~] = freqz(bHQ, aHQ, nfft);
[HaQ,~] = freqz(bHaQ, aHaQ, nfft);

% Magnitude (dB)
HdB = 20*log10(abs(H));
HQdB = 20*log10(abs(HQ));
HaQdB = 20*log10(abs(HaQ));
ymax = ceil(max([HdB; HQdB; HaQdB])) + 2;
figure;
plot(w/pi, HdB,  'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
hold on;
plot(w/pi, HQdB, 'LineWidth', 1.2, 'Color', [0.49 0.18 0.56]);
plot(w/pi, HaQdB,'LineWidth', 1.2, 'Color', [0.22 0.70 0.30]);
grid on;
xlabel('\omega/\pi');
ylabel('Magnitude (dB)');
title('5(f): Magnitude Responses');
ylim([-95 10]);
xlim([0 1.02]);
legend('H','H_Q','H_{aQ}','Location','best');

% Phase
Hph = unwrap(angle(H))*180/pi;
HQph = unwrap(angle(HQ))*180/pi;
HaQph = unwrap(angle(HaQ))*180/pi;
figure;
plot(w/pi, Hph,  'LineWidth', 1.2, 'Color', [0.85 0.33 0.10]);
hold on;
plot(w/pi, HQph, 'LineWidth', 1.2, 'Color', [0.49 0.18 0.56]);
plot(w/pi, HaQph,'LineWidth', 1.2, 'Color', [0.22 0.70 0.30]);
grid on;
xlabel('\omega/\pi');
ylabel('Phase (degrees)');
title('5(f): Phase Responses');
xlim([0 1]);
legend('H','H_Q','H_{aQ}','Location','best');

%% Functions

function mags = plotFkCumulative(sos, g, nfft, ttl, ylims)
if nargin < 5 || isempty(ylims), ylims = [-100 10]; end
L = size(sos,1);
num = g;
den = 1;
mags = zeros(L,1);
C = [0.90 0.10 0.10;
     0.10 0.60 0.80;
     0.60 0.20 0.80;
     0.20 0.75 0.20;
     0.80 0.55 0.00;
     0.20 0.20 0.70];
figure; hold on;
for k = 1:L
    den = conv(den, sos(k,4:6));
    [Hk,w] = freqz(num, den, nfft);
    Habs = abs(Hk);
    mk = max(Habs);
    mags(k) = mk;
    HdB = 20*log10(max(Habs,1e-12)/max(mk,1e-12));
    ck = C(mod(k-1,size(C,1))+1,:);
    plot(w/pi, HdB, 'LineWidth', 1.1, 'Color', ck);
    num = conv(num, sos(k,1:3));
end
yline(0,'k:');
xlim([0 1]); ylim(ylims);
xlabel('\omega/\pi'); ylabel('Magnitude (dB)');
legend(arrayfun(@(i) sprintf('F_%d(z)', i), 1:L, 'uni', 0), 'Location','SouthWest');
title(ttl); hold off;
end

function f = choose_frac_len(v, w)
for f = w-1:-1:0
    lo = -2^(w-1-f);
    hi = 2^(w-1-f) - 2^(-f);
    if all(v >= lo & v <= hi), return, end
end
f = 0;
end

function vq = quantize_twos(v, w, f)
lo = -2^(w-1-f);
hi = 2^(w-1-f) - 2^(-f);
q = round(v * 2^f) / 2^f;
vq = min(max(q, lo), hi);
end
