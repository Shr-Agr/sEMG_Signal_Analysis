clc; clear all; close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reading the data
task_data = readtable('task_ch1Rch2L.csv');
times = task_data.ElapsedTime;
ch1_voltage = task_data.Ch1;
ch2_voltage = task_data.Ch2;

% Plotting the signals without preprocessing
figure(1);
% Subplot 1
subplot(2, 1, 1);
plot(times, ch1_voltage);
title('CH-1 Voltage (without preprocessing)');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
ylim([-1.6e-3 0.1e-3]);
% Subplot 2
subplot(2, 1, 2);
plot(times, ch2_voltage);
title('CH-2 Voltage (without preprocessing)');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
ylim([1.7e-3 3.5e-3]);

% Reading the baseline data
baseline_data = readtable('baseline.csv');
baseline_time = baseline_data.ElapsedTime;
baseline_ch1 = baseline_data.Ch1;
baseline_ch2 = baseline_data.Ch2;

% Plotting the baseline signals without preprocessing
figure(2);
% Subplot 1
subplot(2, 1, 1);
plot(baseline_time, baseline_ch1);
title('CH-1 Voltage Baseline (without preprocessing)');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
%ylim([-1.6e-3 0.1e-3]);
% Subplot 2
subplot(2, 1, 2);
plot(baseline_time, baseline_ch2);
title('CH-2 Voltage Baseline (without preprocessing)');
xlabel('Time (s)');
ylabel('Amplitude (µV)');
%ylim([1.7e-3 3.5e-3]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing
% Step 1: Detrend the EMG signal to remove linear drift
ch1_voltage_detrended = detrend(ch1_voltage(~isnan(ch1_voltage)));
ch2_voltage_detrended = detrend(ch2_voltage(~isnan(ch2_voltage)));
baseline_ch1_voltage_detrended = detrend(baseline_ch1(~isnan(baseline_ch1)));
baseline_ch2_voltage_detrended = detrend(baseline_ch2(~isnan(baseline_ch2)));

% Step 2: Apply a Hanning window
hanning_window_ch1 = hann(length(ch1_voltage_detrended));
hanning_window_ch2 = hann(length(ch2_voltage_detrended));
hanning_window_baseline_ch1 = hann(length(baseline_ch1_voltage_detrended));
hanning_window_baseline_ch2 = hann(length(baseline_ch2_voltage_detrended));
windowed_signal_ch1 = ch1_voltage_detrended .* hanning_window_ch1;
windowed_signal_ch2 = ch2_voltage_detrended .* hanning_window_ch2;
windowed_signal_baseline_ch1 = baseline_ch1_voltage_detrended .* hanning_window_baseline_ch1;
windowed_signal_baseline_ch2 = baseline_ch2_voltage_detrended .* hanning_window_baseline_ch2;

% Step 3: Notch filter to remove supply frequency component of 50 Hz
notch_freq = 50;                                                            % Frequency to remove (e.g., 50 Hz)
notch_Q = 30;                                                               % Quality factor
fs = 1000;
notch_filtered_signal_ch1 = notchFilter(windowed_signal_ch1, fs, notch_freq, notch_Q);
notch_filtered_signal_ch2 = notchFilter(windowed_signal_ch2, fs, notch_freq, notch_Q);
notch_filtered_signal_baseline_ch1 = notchFilter(windowed_signal_baseline_ch1, fs, notch_freq, notch_Q);
notch_filtered_signal_baseline_ch2 = notchFilter(windowed_signal_baseline_ch2, fs, notch_freq, notch_Q);

% Step 4: Bandpass filtering
low_cutoff = 20;               % Low cutoff frequency for bandpass (e.g., 20 Hz)
high_cutoff = 450;             % High cutoff frequency for bandpass (e.g., 450 Hz)
filter_order = 4;              % Filter order for bandpass filter
bandpass_filtered_signal_ch1 = bandpassFilter(notch_filtered_signal_ch1, fs, low_cutoff, high_cutoff, filter_order);
bandpass_filtered_signal_ch2 = bandpassFilter(notch_filtered_signal_ch2, fs, low_cutoff, high_cutoff, filter_order);
bandpass_filtered_signal_baseline_ch1 = bandpassFilter(notch_filtered_signal_baseline_ch1, fs, low_cutoff, high_cutoff, filter_order);
bandpass_filtered_signal_baseline_ch2 = bandpassFilter(notch_filtered_signal_baseline_ch2, fs, low_cutoff, high_cutoff, filter_order);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting the preprocessed signals and doing fft (fft is done after zero padding) for ch1

N1 = length(bandpass_filtered_signal_ch1);
nfft1 = 2^nextpow2(N1);
f1 = (0:nfft1-1)*(fs/nfft1);
fft_ch1=fft(bandpass_filtered_signal_ch1,nfft1);

% Plot the bandpass-filtered signal in the time domain
time_adjusted1 = linspace(min(times), max(times), length(bandpass_filtered_signal_ch1));
figure(3);
subplot(2, 1, 1);
plot(time_adjusted1, bandpass_filtered_signal_ch1);
title('Preprocessed Signal (CH1) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude (µV)');

% Compute the one-sided amplitude spectrum and frequency axis for plotting
one_sided_amplitude_ch1 = abs(fft_ch1(1:nfft1/2+1)) / nfft1;
f_one_sided_ch1 = f1(1:nfft1/2+1);

% Plot the FFT of the bandpass-filtered signal (frequency domain)
subplot(2, 1, 2);
plot(f_one_sided_ch1, one_sided_amplitude_ch1);
title('FFT of Preprocessed Signal (CH1)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


% Plotting the preprocessed signals and doing fft (fft is done after zero padding) for ch2

N2 = length(bandpass_filtered_signal_ch2);
nfft2 = 2^nextpow2(N2);
f2 = (0:nfft2-1)*(fs/nfft2);
fft_ch2=fft(bandpass_filtered_signal_ch2,nfft2);

% Plot the bandpass-filtered signal in the time domain
time_adjusted2 = linspace(min(times), max(times), length(bandpass_filtered_signal_ch2));
figure(4);
subplot(2, 1, 1);
plot(time_adjusted2, bandpass_filtered_signal_ch2);
title('Preprocessed Signal (CH2) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude (µV)');

% Compute the one-sided amplitude spectrum and frequency axis for plotting
one_sided_amplitude_ch2 = abs(fft_ch2(1:nfft2/2+1)) / nfft2;
f_one_sided_ch2 = f2(1:nfft2/2+1);

% Plot the FFT of the bandpass-filtered signal (frequency domain)
subplot(2, 1, 2);
plot(f_one_sided_ch2, one_sided_amplitude_ch2);
title('FFT of Preprocessed Signal (CH2)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


% Plotting the baseline preprocessed signals and doing fft (fft is done after zero padding) for ch1

N3 = length(bandpass_filtered_signal_baseline_ch1);
nfft3 = 2^nextpow2(N3);
f3 = (0:nfft3-1)*(fs/nfft3);
fft_baseline_ch1=fft(bandpass_filtered_signal_baseline_ch1,nfft3);

% Plot the bandpass-filtered signal in the time domain
time_adjusted3 = linspace(min(times), max(times), length(bandpass_filtered_signal_baseline_ch1));
figure(5);
subplot(2, 1, 1);
plot(time_adjusted3, bandpass_filtered_signal_baseline_ch1);
title('Baseline Preprocessed Signal (CH1) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude (µV)');

% Compute the one-sided amplitude spectrum and frequency axis for plotting
one_sided_amplitude_baseline_ch1 = abs(fft_baseline_ch1(1:nfft3/2+1)) / nfft3;
f_one_sided_baseline_ch1 = f3(1:nfft3/2+1);

% Plot the FFT of the bandpass-filtered signal (frequency domain)
subplot(2, 1, 2);
plot(f_one_sided_baseline_ch1, one_sided_amplitude_baseline_ch1);
title('FFT of Baseline Preprocessed Signal (CH1)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');


% Plotting the baseline preprocessed signals and doing fft (fft is done after zero padding) for ch2

N4 = length(bandpass_filtered_signal_baseline_ch2);
nfft4 = 2^nextpow2(N4);
f4 = (0:nfft4-1)*(fs/nfft4);
fft_baseline_ch2=fft(bandpass_filtered_signal_baseline_ch2,nfft4);

% Plot the bandpass-filtered signal in the time domain
time_adjusted4 = linspace(min(times), max(times), length(bandpass_filtered_signal_baseline_ch2));
figure(6);
subplot(2, 1, 1);
plot(time_adjusted4, bandpass_filtered_signal_baseline_ch2);
title('Baseline Preprocessed Signal (CH2) in Time Domain');
xlabel('Time (s)');
ylabel('Amplitude (µV)');

% Compute the one-sided amplitude spectrum and frequency axis for plotting
one_sided_amplitude_baseline_ch2 = abs(fft_baseline_ch2(1:nfft4/2+1)) / nfft4;
f_one_sided_baseline_ch2 = f4(1:nfft4/2+1);

% Plot the FFT of the bandpass-filtered signal (frequency domain)
subplot(2, 1, 2);
plot(f_one_sided_baseline_ch2, one_sided_amplitude_baseline_ch2);
title('FFT of Baseline Preprocessed Signal (CH2)');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notch Filter
function notch_filtered_signal = notchFilter(input_signal, fs, f0, Q)
    % notchFilter applies a notch filter to remove a specific frequency component
    %
    % Parameters:
    % input_signal - the signal to be filtered
    % fs           - sampling frequency of the input signal
    % f0           - frequency to remove (e.g., 50 Hz)
    % Q            - quality factor (determines the notch bandwidth)
    %
    % Returns:
    % filtered_signal - the signal after notch filtering
    
    % Design the Notch Filter
    wo = f0 / (fs / 2);  % Normalize the frequency
    bw = wo / Q;  % Bandwidth of the notch filter
    [b, a] = iirnotch(wo, bw);  % Create notch filter coefficients
    
    % Apply the filter to the input signal
    notch_filtered_signal = filter(b, a, input_signal);
end

% Bandpass Filter
function bandpass_filtered_signal = bandpassFilter(input_signal, fs, low_cutoff, high_cutoff, filter_order)
    % bandpassFilter applies a bandpass filter to a signal
    %
    % Parameters:
    % input_signal  - the signal to be filtered
    % fs            - sampling frequency of the input signal
    % low_cutoff    - low cutoff frequency of the bandpass filter (e.g., 20 Hz)
    % high_cutoff   - high cutoff frequency of the bandpass filter (e.g., 450 Hz)
    % filter_order  - order of the Butterworth filter
    %
    % Returns:
    % filtered_signal - the signal after bandpass filtering
    
    % Design the Bandpass Filter
    [b, a] = butter(filter_order, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
    
    % Apply the filter to the input signal
    bandpass_filtered_signal = filter(b, a, input_signal);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%