
clear all;
close all;
% Parameters (change these as per your setup
task_data_ch2 = readtable('bandpass_filtered_signal_ch2.csv');
bandpass_filtered_signal_ch2=task_data_ch2.bandpass_filtered_signal_ch2;
fs = 1000;  % Sampling frequency in Hz
t = (0:length(bandpass_filtered_signal_ch2)-1) / fs;  % Time vector

% Step 1: Segmenting the MVC Period
% Define MVC start and end times in seconds (adjust based on your data)
mvc_start = 1;   % Start time of MVC in seconds
mvc_end = 41;     % End time of MVC in seconds

% Convert time to indices
mvc_start_idx = round(mvc_start * fs);
mvc_end_idx = round(mvc_end * fs);

% Extract the MVC segment
mvc_segment = bandpass_filtered_signal_ch2(mvc_start_idx:mvc_end_idx);
mvc_time = t(mvc_start_idx:mvc_end_idx);  % Time vector for MVC segment

% Step 2: Full-Wave Rectification
rectified_sEMG = abs(mvc_segment);

% Step 3: Envelope Extraction using Low-Pass Filtering
% Define cutoff frequency for the low-pass filter (e.g., 5-10 Hz)
cutoff_freq = 10;  % Cutoff frequency in Hz

% Design a low-pass filter (Butterworth filter of order 4)
[b, a] = butter(4, cutoff_freq / (fs / 2), 'low');
sEMG_envelope1 = filtfilt(b, a, rectified_sEMG);

% Plotting Results
figure;
subplot(3, 1, 1);
plot(mvc_time, mvc_segment);
title('sEMG Signal Non-Dominant');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3, 1, 2);
plot(mvc_time, rectified_sEMG);
title('Rectified sEMG Signal Non-Dominant');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(3, 1, 3);
plot(mvc_time, sEMG_envelope1);
title('sEMG Envelope (Low-Pass Filtered) Non-Dominant');
xlabel('Time (s)');
ylabel('Amplitude');

% Step 4: Calculate Additional Parameters

% Integrated EMG (IEMG)
iemg = trapz(mvc_time, rectified_sEMG);

% Mean Absolute Value (MAV)
mav = mean((rectified_sEMG));

% Root Mean Square (RMS)
rms_value = sqrt(mean(rectified_sEMG.^2));

% Power Spectral Density (PSD) using Welch's method
[pxx, f] = pwelch(rectified_sEMG, [], [], [], fs);

% Mean Frequency (MNF)
mnf = sum(f .* pxx) / sum(pxx);

% Median Frequency (MDF)
cumulative_power = cumsum(pxx);
mdf_idx = find(cumulative_power >= sum(pxx) / 2, 1);
mdf = f(mdf_idx);


% Step 5: Display Results
fprintf('Integrated EMG (IEMG): %f\n', iemg);
fprintf('Mean Absolute Value (MAV): %f\n', mav);
fprintf('Root Mean Square (RMS): %f\n', rms_value);
fprintf('Mean Frequency (MNF): %f Hz\n', mnf);
fprintf('Median Frequency (MDF): %f Hz\n', mdf);

% Step 6: Plot PSD
figure;
plot(f, 10 * log10(pxx));
title('Power Spectral Density (PSD)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;
% Signal Variance (VAR)
signal_variance = var(rectified_sEMG);

% Zero Crossing Rate (ZCR)
zcr = sum(diff(rectified_sEMG > 0));

% Slope Sign Changes (SSC)
ssc = sum(diff(sign(diff(rectified_sEMG))) ~= 0);

% Waveform Length (WL)
waveform_length = sum(abs(diff(rectified_sEMG)));
% Display Results
fprintf('Signal Variance (VAR): %f\n', signal_variance*1e6);
fprintf('Zero Crossing Rate (ZCR): %d\n', zcr);
fprintf('Slope Sign Changes (SSC): %d\n', ssc);
fprintf('Waveform Length (WL): %f\n', waveform_length);
