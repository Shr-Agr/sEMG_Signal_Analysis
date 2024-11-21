clc;

% Simulate a sample signal (replace with actual signal)
fs = 1000;            % Sampling frequency in Hz
t = 0:1/fs:1-1/fs;    % Time vector
signal = bandpass_filtered_signal_ch2; % Replace with actual signal

% Parameters
segment_length = 512; % Length of each segment for FFT

% Frequency Analysis
nfft = segment_length;                % Number of FFT points
frequencies = (0:nfft/2-1)*(fs/nfft); % Frequency axis (one-sided)
segments = buffer(signal, segment_length, 0, 'nodelay'); % Divide into segments
num_segments = size(segments, 2);

% Initialize metrics
mpf = zeros(1, num_segments); % Median Power Frequency
mnf = zeros(1, num_segments); % Mean Power Frequency
mav = zeros(1, num_segments); % Mean Absolute Value 
rms = zeros(1, num_segments); % Root Mean Square

for i = 1:num_segments
    segment = segments(:, i); % Current segment
    segment = segment - mean(segment); % Detrend (remove DC component)
    
    % MAV and RMS calculations
    mav(i) = mean(abs(segment)); % Mean Absolute Value
    rms(i) = sqrt(mean(segment.^2)); % Root Mean Square
    
    % FFT and PSD
    fft_segment = fft(segment, nfft);
    psd = (abs(fft_segment(1:nfft/2)).^2) / (fs * nfft); % One-sided PSD
    
    % Ensure matching dimensions
    frequencies = frequencies(:); % Make frequencies a column vector
    psd = psd(:);                 % Make PSD a column vector
    
    % MPF Calculation
    cumulative_power = cumsum(psd); % Cumulative sum of PSD
    total_power = cumulative_power(end);
    mpf(i) = frequencies(find(cumulative_power >= total_power/2, 1)); % Median Power Frequency
    
    % MNF Calculation
    mnf(i) = sum(frequencies .* psd) / sum(psd); % Mean Power Frequency
end

% Display results
disp('Median Power Frequency (MPF) for each segment:');
disp(mpf);
disp('Mean Power Frequency (MNF) for each segment:');
disp(mnf);
disp('Mean Absolute Value (MAV) for each segment:'); 
disp(mav); 
disp('Root Mean Square (RMS) for each segment:'); 
disp(rms);

% Plotting PSD of the first segment
figure;
plot(frequencies, psd);
title('Power Spectral Density (First Segment)');
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
grid on;

% Plot MNF, MPF, MAV, and RMS across segments
figure;

% Plot MNF
subplot(4, 1, 1);
plot(1:num_segments, mnf, '-o');
title('Mean Power Frequency (MNF) Across Segments');
xlabel('Segment Number');
ylabel('MNF (Hz)');
grid on;

% Plot MPF
subplot(4, 1, 2);
plot(1:num_segments, mpf, '-o');
title('Median Power Frequency (MPF) Across Segments');
xlabel('Segment Number');
ylabel('MPF (Hz)');
grid on;

% Plot MAV
subplot(4, 1, 3);
plot(1:num_segments, mav, '-o');
title('Mean Absolute Value (MAV) Across Segments');
xlabel('Segment Number');
ylabel('MAV');
grid on;

% Plot RMS
subplot(4, 1, 4);
plot(1:num_segments, rms, '-o');
title('Root Mean Square (RMS) Across Segments');
xlabel('Segment Number');
ylabel('RMS');
grid on;