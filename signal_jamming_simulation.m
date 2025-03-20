%% Generate DT LPE signal %%
rf_freq = 1.0075e9:5e6:1.0925e9; % Frequencies in Hz
fq_select = rf_freq(randi(length(rf_freq)));   %%% randomly select one of the freq
t=0:1:8191; % create 8192 DT sample index points 
% for jam scenario 1 samples can be low but with to detect signals with
% white noise the samples were increased 
N = length(t);
%%%%
% we have our our 1.0GHz to 1.1GHz shifted to -pi to pi
% so all our freq will get shifted accordingly
% scaling factor is DT bandwidth (2π) divided by the RF bandwidth (0.1 GHz)
% to get equivalent DT frequency for a given RF frequency
% we subtract the base frequency (1.0 GHz) from it and then scale
% basically each freq becomes ---> -pi + freq values
% so if we have 1.0075 and we know 1.0 maps to -pi 
% then similarly 1.0075 gets to -pi + 1.0075(just scaled)
%%%%%
omega = -pi + (2*pi / 0.1e9) * (fq_select - 1.0e9);
fprintf('RF Frequency: %.4f GHz, DT Frequency: %.4f rad/sample\n', fq_select/1e9, omega);
% Create a complex signal
I = cos(omega * t); % in-phase component
Q = 1i * sin(omega * t); % quadrature component
x_l = I + Q;

%{
%PLOT TO CHECK IF WE SEE THE OUR SIGNAL BETWEEN -pi to pi as we want
N_zp = 10 * 4095;
X = fftshift(fft(x_l, N_zp));
Omega = (-(N_zp/2):((N_zp/2)-1)) * 2 * pi / N_zp; % DT Frequency -pi to pi
h=plot(Omega,20*log10(abs(X)));
set(h,'linewidth',2)
set(gca,'fontsize',12)
xlabel('Omega rad/sample')
ylabel('|X[k]|')
title('DT LPE Spectrum(from -pi to pi)')
grid
xlim([-pi, pi]);
%}

%% JAMMING Scenario I: Sum of Sinusoids %%

jam_freqs = 1.005e9:5e6:1.095e9; % Jamming frequencies in Hz
x_jamming = zeros(1, N); % Initialize jamming signal

for f = jam_freqs
    omega_j = -pi + (2*pi / 0.1e9) * (f - 1.0e9); % Convert to DT frequency
    %fprintf('Jamming Frequency: %.4f GHz, DT Frequency: %.2f rad/sample\n', f/1e9, omega);
    x_jamming = x_jamming + cos(omega_j * t) + 1i * sin(omega_j * t);
end

x_jamming = x_jamming * sqrt(10^(60/10)); % jamming signal is 60 dB higher in power than the desired signal; amplified by a factor of 1000
%{
%PLOT TO WE SEE THE OUR SIGNAL
N_zp = 10 * 4095;
X = fftshift(fft(x_jamming, N_zp));
Omega = (-(N_zp/2):((N_zp/2)-1)) * 2 * pi / N_zp; % DT Frequency -pi to pi
h=plot(Omega,20*log10(abs(X)));
set(h,'linewidth',2)
set(gca,'fontsize',12)
xlabel('Omega rad/sample')
ylabel('|X[k]|')
title('Jamming Signal(Sum of 19 Sinusoids) Spectrum')
grid
xlim([-pi, pi]);
%}

%% JAMMING Scenario II: White Noise Jamming %%

%%%%%%%
%%% SNR in dB is 10log_10 (P_signal/P_noise)
%%% SNR of -30 dB means that the noise is 30 dB stronger than the signal
%%% Since P_noise is 30dB higher than P_signal
%%% The variance of noise will be 10^(30/100) = 1000 times the power of desired signal!
%%%%%%%
noise_v = 1000; 
% Generate complex Gaussian white noise
white_noise = sqrt(noise_v/2) * randn(1, N) + 1i * sqrt(noise_v/2) * randn(1, N);

%{
%plot noise
N_zp = 10 * 4096;
X = fftshift(fft(white_noise, N_zp));
Omega = (-(N_zp/2):((N_zp/2)-1)) * 2 * pi / N_zp; % DT Frequency -pi to pi
h=plot(Omega,20*log10(abs(X)));
set(h,'linewidth',2)
set(gca,'fontsize',12)
xlabel('Omega rad/sample')
ylabel('|X[k]|')
title('DT LPE Spectrum(from -pi to pi)')
grid
xlim([-pi, pi]);
%}

%% DFT ANALYSIS FOR PART I %%

x_l_jammed = x_l + x_jamming;
N_zp_j = 2 * N ;
Omega = (-(N_zp_j/2):((N_zp_j/2)-1)) * 2 * pi / N_zp_j; % DT Frequency -pi to pi
% Different windows
window_types = {'han', 'ham', 'black'};
window_functions = {hann(length(x_l_jammed)), hamming(length(x_l_jammed)), blackman(length(x_l_jammed))};
% Initialize structure for storing FFT results
X_windows = struct();

for i = 1:length(window_functions)
    % Apply window
    x_windowed = x_l_jammed .* window_functions{i}';
    X_windowed = fftshift(fft(x_windowed, N_zp_j));
    % Store in structure
    X_windows.(window_types{i}) = X_windowed;
end


% Plotting
h = plot(Omega,20*log10(abs(X_windows.han)));
set(h,'linewidth',2)
set(gca,'fontsize',12)
xlabel('Omega (rad/sample)')
ylabel('|X[k]| (dB)')
title('Hann Window Jamming Scenario I')
grid
xlim([-pi, pi]);
[pks, locs] = findpeaks(20*log10(abs(X_windows.han)), 'MinPeakProminence',72); % Find Peaks
[min_peak, min_idx] = min(pks);
min_peak_loc = Omega(locs(min_idx)); % Identifying the lowest peak
hold on;
plot(min_peak_loc, min_peak, 'p', 'MarkerSize', 15, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'yellow');
hold off;
legend('DFT Spectrum', 'True Signal', 'Location', 'best');


%% DFT ANALYSIS FOR PART II %%

x_l_noisy = x_l + white_noise;

% averaging DFTs
segment_length = 4096; % Length of each segment
overlap = 2048; % 50% overlap
N_zp = N; 
start_idx = 1;
average_spectrum_dft = zeros(1, N_zp);
num_segments = 0;

% Averaging DFTs with hamming window
while start_idx + segment_length - 1 <= N
    segment = x_l_noisy(start_idx:start_idx + segment_length - 1);
    windowed_segment = segment .* hamming(length(segment))';
    segment_dft = fftshift(fft(windowed_segment, N_zp));
    average_spectrum_dft = average_spectrum_dft + abs(segment_dft);
    start_idx = start_idx + segment_length - overlap;
    num_segments = num_segments + 1;
end

average_spectrum_dft = average_spectrum_dft / num_segments;
Omega_n = (-(N_zp/2):((N_zp/2)-1)) * 2 * pi / N_zp; % DT Frequency -pi to pi
% Plot the DFT spectrum
plot(Omega_n, 20*log10(abs(average_spectrum_dft)), 'b');
set(gca, 'linewidth', 2, 'fontsize', 12);
xlabel('Omega (rad/sample)');
ylabel('|X[k]| (dB)');
title('Average DFT of Noisy Signal');
grid on;
xlim([-pi, pi]);
hold on;
%%% just to visualize where actually our signal is I used real omega value
% omega is unknown to user... for real situation I try to find the highest peak 
% and calculated a percentage of times our signal is actually the highest peak... 
% 33% times our signal has highest peak...code is shown in comment below..
[~, signal_bin] = min(abs(Omega_n - omega));
peak_value = 20*log10(abs(average_spectrum_dft(signal_bin)));
plot(Omega_n(signal_bin), peak_value, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
legend('DFT Spectrum', 'True Signal', 'Location', 'best');
hold off;

%{
% Output SNR
[~, signal_bin] = min(abs(Omega_n - omega)); 
signal_power = abs(average_spectrum_dft(signal_bin))^2;
bins = 1:length(Omega_n); % Array of all bin indices
is_noise_bin = true(size(bins));
is_noise_bin(signal_bin) = false;
noise_bins = Omega_n(is_noise_bin);
noise_power = mean(abs(average_spectrum_dft(noise_bins)).^2);
output_SNR_db = 10 * log10(signal_power / noise_power);
fprintf('Output SNR: %f dB\n', output_SNR_db);


num_trials = 10000; % Number of trials
num_true_positives = 0;
num_false_alarms = 0;
threshold_dB = -30; %noise floor assuming

for trial = 1:num_trials
    rf_freq = 1.0075e9:5e6:1.0925e9; 
    t=0:1:8191; 
    N = length(t);
    f = rf_freq(randi(length(rf_freq)));
    omega = -pi + (2*pi / 0.1e9) * (f - 1.0e9);
    x_l = cos(omega * t) + 1i * sin(omega * t);

    noise_v = 1000; 
    r_noise = sqrt(noise_v/2) * randn(1, N);
    i_noise = sqrt(noise_v/2) * randn(1, N);
    white_noise = r_noise + 1i * i_noise;

    x_l_noisy = x_l + white_noise;

    segment_length = 4096; 
    overlap = 2048;
    N_zp = 7 * N; 
    start_idx = 1;
    average_spectrum_dft = zeros(1, N_zp);
    num_segments = 0;

    while start_idx + segment_length - 1 <= N
        segment = x_l_noisy(start_idx:start_idx + segment_length - 1);
        windowed_segment = segment .* hamming(length(segment))';
        segment_dft = fftshift(fft(windowed_segment, N_zp));
        average_spectrum_dft = average_spectrum_dft + abs(segment_dft);
        start_idx = start_idx + segment_length - overlap;
        num_segments = num_segments + 1;
    end

    average_spectrum_dft = average_spectrum_dft / num_segments;
    threshold = 1000; % Threshold for -30dB
    thresholded_spectrum = average_spectrum_dft;
    thresholded_spectrum(average_spectrum_dft < threshold) = 0;

    Omega = (-(N_zp/2):((N_zp/2)-1)) * 2 * pi / N_zp; % DT Frequency -pi to pi


    % Find peaks in the DFT output
    peaks = findpeaks(20*log10(abs(average_spectrum_dft)));
    if isempty(peaks)
        continue; % No peaks detected, skip to next iteration
    end

    [max_peak, idx] = max(peaks);
    peak_frequency = Omega(idx);
    if max_peak > threshold_dB
        if abs(peak_frequency - omega) < 2 % Define a tolerance
            num_true_positives = num_true_positives + 1;
        else
            num_false_alarms = num_false_alarms + 1;
        end
    end
end

Pd = num_true_positives / num_trials;
FAR = num_false_alarms / num_trials;
fprintf('Detection Probability: %.2f\n', Pd);
fprintf('False Alarm Rate: %.2f\n', FAR);
%}
