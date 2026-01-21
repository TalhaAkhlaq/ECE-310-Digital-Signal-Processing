%% Talha Akhlaq
% ECE310 Digital Signal Processing
% Problem Set 4
clear; clc; close all;

%% Problem 6(a)
fprintf('Part (A): \n');
w_name = 'db1';

[h_L, h_H, f_L, f_H] = wfilters(w_name);

fprintf('Analysis Low H0:    %s\n', mat2str(h_L, 4));
fprintf('Analysis High H1:   %s\n', mat2str(h_H, 4));
fprintf('Synthesis Low F0:   %s\n', mat2str(f_L, 4));
fprintf('Synthesis High F1:  %s\n', mat2str(f_H, 4));

pr_res = conv(h_L, f_L) + conv(h_H, f_H);
fprintf('Convolution Result: %s\n', mat2str(pr_res, 4));

if abs(pr_res(2) - 2) < 1e-10
    fprintf('Delay=1, Scale=2.\n\n');
end

%% Problem 6(b)
fprintf('Part (B): \n');
orders = [4, 8];
pts = 10000;
w = linspace(0, pi, pts);

for k = 1:length(orders)
    N = orders(k);
    cur_wav = sprintf('db%d', N);
    fprintf('Wavelet: %s (N=%d) \n', cur_wav, N);

    %1
    [h0, h1, f0, f1] = wfilters(cur_wav);

    %2
    H0_mag = freqz(h0, 1, w);
    H1_mag = freqz(h1, 1, w);

    figure('Name', [cur_wav ' Frequency Response']);
    plot(w/pi, abs(H0_mag), 'b', 'LineWidth', 1.5); hold on;
    plot(w/pi, abs(H1_mag), 'r--', 'LineWidth', 1.5);
    title(['Magnitude Response: ', cur_wav]);
    xlabel('Normalized Frequency \omega (rad)'); 
    ylabel('Magnitude');
    legend('|H_0|', '|H_1|'); grid on; hold off;

    %3
    pow_sum = abs(H0_mag).^2 + abs(H1_mag).^2;
    err_pwr = max(abs(pow_sum - 2));
    fprintf('3. %.2e\n', err_pwr);

    %4
    E00 = h0(1:2:end); E01 = h0(2:2:end);
    E10 = h1(1:2:end); E11 = h1(2:2:end);

    R00 = f0(1:2:end); R01 = f0(2:2:end);
    R10 = f1(1:2:end); R11 = f1(2:2:end);

    %5
    E00_rev = fliplr(E00);
    idx_nz = find(abs(E00_rev) > 1e-6, 1);
    if ~isempty(idx_nz)
        c_val = R00(idx_nz) / E00_rev(idx_nz);
    else
        c_val = 0;
    end
    fprintf('5. %.4f\n', c_val);

    %6
    P00 = conv(R00, E00) + conv(R01, E10);
    P01 = conv(R00, E01) + conv(R01, E11);
    P10 = conv(R10, E00) + conv(R11, E10);

    err_off = max([max(abs(P01)), max(abs(P10))]);
    [diag_val, diag_idx] = max(abs(P00));
    sys_M = diag_idx - 1;

    fprintf('6. Off-Diagonal Error: %.2e\n', err_off);
    fprintf('   Diagonal Peak (d): %.4f\n', P00(diag_idx));
    fprintf('   Delay M: %d\n', sys_M);

    %7
    n_seq = 0:(length(h0)-1);
    alt_v = (-1).^n_seq;
    
    T_sys = 0.5 * (conv(f0, h0) + conv(f1, h1));
    A_sys = 0.5 * (conv(f0, h0.*alt_v) + conv(f1, h1.*alt_v));

    [T_amp, T_pos] = max(abs(T_sys));
    fprintf('7. T(z) Peak: %.4f (location %d)\n', T_amp, T_pos);
    fprintf('   Max Aliasing A(z): %.2e\n', max(abs(A_sys)));

    %8
    H0_sq = abs(H0_mag).^2;
    dw = w(2) - w(1);

    d1 = diff(H0_sq) / dw; 
    d1 = [d1, 0];
    d2 = diff(H0_sq, 2) / (dw^2); 
    d2 = [d2, 0, 0];

    fprintf('8. 1st Derivative [0, pi]: [%.2e, %.2e]\n', d1(1), d1(end));
    fprintf('   2nd Derivative [0, pi]: [%.2e, %.2e]\n', d2(1), d2(end));

    figure('Name', [cur_wav ' Derivatives']);
    subplot(2,1,1); plot(w/pi, d1, 'b', 'LineWidth', 1.2); grid on;
    title(['First Derivative of |H_0|^2: ', cur_wav]); 
    ylabel('First Derivative');
    xlabel('Normalized Frequency \omega (rad)');

    subplot(2,1,2); plot(w/pi, d2, 'm', 'LineWidth', 1.2); grid on;
    title(['Second Derivative of |H_0|^2: ', cur_wav]); 
    ylabel('Second Derivative');
    xlabel('Normalized Frequency \omega (rad)');

    %9
    h0_u2 = zeros(1, 2*length(h0)); h0_u2(1:2:end) = h0;
    h1_u2 = zeros(1, 2*length(h1)); h1_u2(1:2:end) = h1;

    G_L2 = cell(4,1);
    G_L2{1} = conv(h0, h0_u2);
    G_L2{2} = conv(h0, h1_u2);
    G_L2{3} = conv(h1, h0_u2);
    G_L2{4} = conv(h1, h1_u2);

    figure('Name', [cur_wav ' 2-Stage']);
    hold on;
    clrs = lines(4);
    for i = 1:4
        [G_out, wg] = freqz(G_L2{i}, 1, pts);
        plot(wg/pi, abs(G_out), 'Color', clrs(i,:), 'LineWidth', 1.5);
    end
    title(['2-Stage Tree: ', cur_wav]); 
    ylabel('Magnitude');
    xlabel('Normalized Frequency \omega (rad)');
    grid on; hold off; 

    h0_u4 = zeros(1, 4*length(h0)); h0_u4(1:4:end) = h0;
    h1_u4 = zeros(1, 4*length(h1)); h1_u4(1:4:end) = h1;
    G_L3 = cell(8,1);

    for i = 1:4
        G_L3{2*i-1} = conv(G_L2{i}, h0_u4);
        G_L3{2*i}   = conv(G_L2{i}, h1_u4);
    end

    figure('Name', [cur_wav ' 3-Stage']);
    hold on;
    clrs8 = hsv(8);
    for i = 1:8
        [G3_out, wg3] = freqz(G_L3{i}, 1, pts);
        plot(wg3/pi, abs(G3_out), 'Color', clrs8(i,:), 'LineWidth', 1.2);
    end
    title(['3-Stage Tree: ', cur_wav]); 
    ylabel('Magnitude');
    xlabel('Normalized Frequency \omega (rad)');
    grid on; hold off;

    fprintf('\n');
end