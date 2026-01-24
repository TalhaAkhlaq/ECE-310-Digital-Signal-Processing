%% Talha Akhlaq
% ECE310 Digital Signal Processing
% Problem Set 5
clear; clc; close all;

%% Problem 2(a)
h_Lap = (1/6) * [1 4 1; 4 -20 4; 1 4 1];

hx = (1/4) * [ 1  0 -1;
              2  0 -2;
              1  0 -1];

hy = (1/4) * [-1 -2 -1;
               0  0  0;
               1  2  1];

h_LapSob = conv2(hx, hx) + conv2(hy, hy);

fprintf('Problem 2(a):\n');
disp('h_LapSob coefficients:');
disp(h_LapSob);

sz = size(h_LapSob);
fprintf('Size of h_LapSob: %d x %d\n', sz(1), sz(2));

if isequal(h_LapSob, rot90(h_LapSob, 2))
    disp('h_LapSob is zero-phase (symmetric)');
else
    disp('h_LapSob is not zero-phase');
end

%% Problem 2(b)
[H_Lap, fx_lap, fy_lap] = freqz2(h_Lap);
[H_LapSob, fx_sob, fy_sob] = freqz2(h_LapSob);

H_Lap = real(H_Lap);
H_LapSob = real(H_LapSob);

max_Lap = max(H_Lap(:));
max_LapSob = max(H_LapSob(:));

fprintf('\nProblem 2(b):\n');
fprintf('Maximum value of H_Lap: %.6f\n', max_Lap);
fprintf('Maximum value of H_LapSob: %.6f\n', max_LapSob);

%% Problem 2(c)

figure('Name', 'Problem 2(c): H_Lap');
subplot(1,2,1);
surf(fx_lap, fy_lap, H_Lap);
colormap(jet);
title('Surface Plot of H_{Lap}');
xlabel('f_x'); ylabel('f_y');
shading interp;

subplot(1,2,2);
contour(fx_lap, fy_lap, H_Lap);
colormap(jet);
title('Contour Plot of H_{Lap}');
xlabel('f_x'); ylabel('f_y');
grid on;

figure('Name', 'Problem 2(c): H_LapSob');
subplot(1,2,1);
surf(fx_sob, fy_sob, H_LapSob);
colormap(jet);
title('Surface Plot of H_{LapSob}');
xlabel('f_x'); ylabel('f_y');
shading interp;

subplot(1,2,2);
contour(fx_sob, fy_sob, H_LapSob);
colormap(jet);
title('Contour Plot of H_{LapSob}');
xlabel('f_x'); ylabel('f_y');
grid on;

fprintf('\n');
disp('Problem 2(c):');
disp('HLap: high-pass');
disp('HLapSob: high-pass');
disp('Contours are roughly circular near the origin, so approximately isotropic');

%% Problem 2(d)

lily_data = load('LilyImg.mat');
rodan_data = load('RodanImg.mat');

names_lily = fieldnames(lily_data);
Lily = lily_data.(names_lily{1});

names_rodan = fieldnames(rodan_data);
Rodan = rodan_data.(names_rodan{1});

Lily_Lap = filter2(h_Lap, Lily, 'same');
Lily_LapSob = filter2(h_LapSob, Lily, 'same');

Rodan_Lap = filter2(h_Lap, Rodan, 'same');
Rodan_LapSob = filter2(h_LapSob, Rodan, 'same');

figure('Name', 'Problem 2(d): Lily');
subplot(1,3,1);
image(Lily, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Original Lily');

subplot(1,3,2);
image(Lily_Lap, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Lily with h_{Lap}');

subplot(1,3,3);
image(Lily_LapSob, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Lily with h_{LapSob}');

figure('Name', 'Problem 2(d): Rodan');
subplot(1,3,1);
image(Rodan, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Original Rodan');

subplot(1,3,2);
image(Rodan_Lap, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Rodan with h_{Lap}');

subplot(1,3,3);
image(Rodan_LapSob, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Rodan with h_{LapSob}');

%% Problem 3(b)

Lily_Up = upsample_matrix(double(Lily));
Rodan_Up = upsample_matrix(double(Rodan));

figure('Name', 'Problem 3(b): Upsampled Images');
subplot(1,2,1);
image(Lily_Up, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Upsampled Lily');

subplot(1,2,2);
image(Rodan_Up, 'CDataMapping', 'scaled'); colormap(gray); axis image;
title('Upsampled Rodan');

fprintf('\n');
disp('Problem 3(b):');
disp('10x10 Check:');
disp('Top-left 10x10 Original Lily:');
disp(Lily(1:10, 1:10));
disp('Top-left 10x10 Upsampled Lily:');
disp(Lily_Up(1:10, 1:10));

disp('Top-left 10x10 Original Rodan:');
disp(Rodan(1:10, 1:10));
disp('Top-left 10x10 Upsampled Rodan:');
disp(Rodan_Up(1:10, 1:10));


%% Problem 3(c)

F_Lily = fftshift(fft2(double(Lily)));
F_Lily_Up = fftshift(fft2(Lily_Up));

F_Rodan = fftshift(fft2(double(Rodan)));
F_Rodan_Up = fftshift(fft2(Rodan_Up));

S_Lily = log(1 + abs(F_Lily));
S_Lily_Up = log(1 + abs(F_Lily_Up));
S_Rodan = log(1 + abs(F_Rodan));
S_Rodan_Up = log(1 + abs(F_Rodan_Up));

figure('Name', 'Problem 3(c)');
subplot(2,2,1);
image(S_Lily, 'CDataMapping', 'scaled'); colormap(gray); axis square;
title('Spectrum: Original Lily');

subplot(2,2,2);
image(S_Lily_Up, 'CDataMapping', 'scaled'); colormap(gray); axis square;
title('Spectrum: Upsampled Lily');

subplot(2,2,3);
image(S_Rodan, 'CDataMapping', 'scaled'); colormap(gray); axis square;
title('Spectrum: Original Rodan');

subplot(2,2,4);
image(S_Rodan_Up, 'CDataMapping', 'scaled'); colormap(gray); axis square;
title('Spectrum: Upsampled Rodan');

fprintf('\n');
disp('Problem 3(c):');
disp('Note: Spectra shown are log magnitude');
disp('Imaging Distortion (Spectral Imaging):');
disp('Replicas of the original baseband spectrum appear in high-frequency regions');
disp('The original spectrum is compressed and repeated in the high-frequency regions of the larger frequency domain');



%% Problem 3(a)

function v = upsample_matrix(u)
    [rows, cols] = size(u);
    v = zeros(2*rows, 2*cols);
    v(1:2:end, 1:2:end) = u;
end
