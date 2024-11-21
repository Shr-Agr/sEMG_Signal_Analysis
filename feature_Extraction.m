% Define sampling frequency
fs = 1000;  % Hz

% Assuming you have these variables from previous steps:
% sEMG_dominant - rectified & filtered envelope for dominant hand
% sEMG_nondominant - rectified & filtered envelope for non-dominant hand
task_data_ch1 = readtable('ch1_MVC_envelope.csv');
sEMG_dominant=task_data_ch1.ch1_MVC_envelope;

task_data_ch2=readtable('ch2_MVC_envelope.csv');
sEMG_nondominant=task_data_ch2.ch2_MVC_envelope;

%% Step 4: Feature Extraction

% Time-Domain Features
MAV_dominant = mean(sEMG_dominant);    % Mean Absolute Value (MAV)
MAV_nondominant = mean(sEMG_nondominant);

RMS_dominant = rms(sEMG_dominant);     % Root Mean Square (RMS)
RMS_nondominant = rms(sEMG_nondominant);

IEMG_dominant = trapz(sEMG_dominant);  % Integrated EMG (IEMG)
IEMG_nondominant = trapz(sEMG_nondominant);

% Frequency-Domain Features
% Use FFT for frequency analysis
N = length(sEMG_dominant);  % Number of points
nfft1=2^nextpow2(N);
f = (0:nfft1-1)*(fs/nfft1);     % Frequency vector

% Dominant hand FFT
Y_dominant = fft(sEMG_dominant,nfft1);
P_dominant = abs(Y_dominant).^2/(nfft1*fs);  % Power spectral density
P_dominant = P_dominant(1:nfft1/2+1);   % Take half of the symmetric spectrum
f_half1 = f(1:nfft1/2+1);                % Corresponding frequencies

% Mean and Median Frequency for Dominant
MNF_dominant = sum(f_half1 .* P_dominant) / sum(P_dominant);  % Mean Frequency
MDF_dominant = f_half1(find(cumsum(P_dominant) >= 0.5 * sum(P_dominant), 1));  % Median Frequency

% Non-dominant hand FFT
% Use FFT for frequency analysis
N2 = length(sEMG_nondominant);  % Number of points
nfft2=2^nextpow2(N2);
f = (0:nfft2-1)*(fs/nfft2);     % Frequency vector

Y_nondominant = fft(sEMG_nondominant,nfft2);
P_nondominant = abs(Y_nondominant).^2/(nfft2*fs);
P_nondominant = P_nondominant(1:nfft2/2+1);
f_half2 = f(1:nfft2/2+1);
% Mean and Median Frequency for Non-Dominant
MNF_nondominant = sum(f_half2 .* P_nondominant) / sum(P_nondominant);
MDF_nondominant = f_half2(find(cumsum(P_nondominant) >= 0.5 * sum(P_nondominant), 1));

% Dominant Hand Results
disp(['Dominant Hand - MAV: ', num2str(MAV_dominant, '%.3f'), ...
      ', RMS: ', num2str(RMS_dominant, '%.3f'), ...
      ', IEMG: ', num2str(IEMG_dominant, '%.3f'), ...
      ', MNF: ', num2str(MNF_dominant, '%.3f'), ' Hz', ...
      ', MDF: ', num2str(MDF_dominant, '%.3f'), ' Hz']);

% Non-Dominant Hand Results
disp(['Non-Dominant Hand - MAV: ', num2str(MAV_nondominant, '%.3f'), ...
      ', RMS: ', num2str(RMS_nondominant, '%.3f'), ...
      ', IEMG: ', num2str(IEMG_nondominant, '%.3f'), ...
      ', MNF: ', num2str(MNF_nondominant, '%.3f'), ' Hz', ...
      ', MDF: ', num2str(MDF_nondominant, '%.3f'), ' Hz']);

%% Step 5: Statistical Analysis

% Two-sample t-test for comparing features between dominant and non-dominant hands
% For example, comparing RMS values
[h, p] = ttest2(sEMG_dominant, sEMG_nondominant);
if h == 1
    fprintf('Significant difference found in RMS between dominant and non-dominant hand (p = %.3f).\n', p);
else
    fprintf('No significant difference found in RMS between dominant and non-dominant hand (p = %.3f).\n', p);
end

%% Step 6: Muscle Fatigue Analysis

% Muscle Fatigue Analysis over Time
% Split the signal into smaller windows and calculate the median frequency for each window
window_size = fs;  % Define window size (1-second windows in this case)
num_windows = floor(N / window_size);

% Initialize arrays to store median frequencies over time
MDF_dominant_windows = zeros(1, num_windows);
MDF_nondominant_windows = zeros(1, num_windows);

for k = 1:num_windows
    % Define window range
    window_start = (k-1) * window_size + 1;
    window_end = k * window_size;
    
    % Dominant hand MDF over time
    Y_window = fft(sEMG_dominant(window_start:window_end));
    P_window = abs(Y_window/window_size).^2;
    P_window = P_window(1:window_size/2+1);
    MDF_dominant_windows(k) = f_half(find(cumsum(P_window) >= 0.5 * sum(P_window), 1));
    
    % Non-Dominant hand MDF over time
    Y_window = fft(sEMG_nondominant(window_start:window_end));
    P_window = abs(Y_window/window_size).^2;
    P_window = P_window(1:window_size/2+1);
    MDF_nondominant_windows(k) = f_half(find(cumsum(P_window) >= 0.5 * sum(P_window), 1));
 
end

% Plotting the median frequency over time for both hands to visualize fatigue
time_vector = (1:num_windows) * (window_size / fs);

figure;
plot(time_vector, MDF_dominant_windows, '-o', 'DisplayName', 'Dominant Hand');
hold on;
plot(time_vector, MDF_nondominant_windows, '-o', 'DisplayName', 'Non-Dominant Hand');
title('Muscle Fatigue Analysis: Median Frequency Over Time');
xlabel('Time (s)');
ylabel('Median Frequency (Hz)');
legend;
grid on;
