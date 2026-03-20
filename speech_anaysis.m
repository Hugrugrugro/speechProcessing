[signal, Fs] = audioread('modulator22.wav'); % replace filename if needed

% Use provided frequencySpectrum function to compute/display power and get duration
% Call with padding to next power of two to use FFT algorithm
[power_fs, duration_fs] = frequencySpectrum(signal, Fs, true);

% Save duration info for later display (used below)
fft_mean = duration_fs;
fft_std = 0;

% Modify sampling frequency information (half the original) and save
Fs_new = Fs/2;
audiowrite('modulator22_halfFs.wav', signal, round(Fs_new));

% Temporal plot of signal amplitude
t = (0:length(signal)-1).' / Fs; % time vector using original Fs
figure;
plot(t, signal);
xlabel('Time (s)');
ylabel('Amplitude');
title('Temporal variations of speech signal');
grid on;

% Compute power spectrum using frequencySpectrum (correct 3-argument call)
window = length(signal);
nfft = max(2^nextpow2(window), window);

% FIX: call frequencySpectrum with its correct signature (signal, fs, pad)
% pad=false here so no zero-padding is applied (DFT of full signal)
[power_dft, ~] = frequencySpectrum(signal, Fs, false);
title('Power Spectrum (DFT via frequencySpectrum)');

% Compute power spectrum using FFT (manual) and measure time
ntrials = 5;
fftTimes = zeros(ntrials,1);
for k = 1:ntrials
    tic;
    X = fft(signal, nfft);
    fftTimes(k) = toc;
end
fft_mean = mean(fftTimes);
fft_std = std(fftTimes);

% If signal is multi-channel, combine channels by averaging
if size(X,2) > 1
    X = mean(X,2);
end

% Compute single-sided power spectral density (PSD) in dB/Hz from FFT result
P2 = (abs(X)/window).^2;          % two-sided power (scaled)
P1 = P2(1:floor(nfft/2)+1);
if nfft > 1
    P1(2:end-1) = 2*P1(2:end-1);  % single-sided
end
f = (0:floor(nfft/2)).' * (Fs / nfft);
figure;
plot(f, 10*log10(P1 + eps));
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
title(sprintf('Power Spectrum (FFT) — mean time: %.6f s, std: %.6f s', fft_mean, fft_std));
grid on;

% Compute and display spectrogram using window_size = 30 ms and step_size = 5 ms
window_size_ms = 30;
step_size_ms = 5;
window_size_samples = round(window_size_ms/1000 * Fs);
step_size_samples = round(step_size_ms/1000 * Fs);
noverlap = window_size_samples - step_size_samples;
nfft_spec = max(256, 2^nextpow2(window_size_samples));

figure;
spectrogram(signal, hamming(window_size_samples), noverlap, nfft_spec, Fs, 'yaxis');
title(sprintf('Spectrogram (window = %d ms, step = %d ms)', window_size_ms, step_size_ms));
colorbar;

% ------------------ New processing per assignment ------------------

% Target sampling frequency
Fs_target = 4000;
% Compute integer downsampling factor (assume Fs is integer multiple or approximate)
M = round(Fs / Fs_target);
Fs_ds = Fs / M;

% 1) Downsample using downsample()
y_downsampled = downsample(signal, M);
audiowrite('modulator22_downsampled_downsample.wav', y_downsampled, round(Fs_ds));

% 2) Downsample using decimate()
% decimate applies an anti-aliasing IIR filter by default
if size(signal,2) > 1
    y_decimated = zeros(ceil(size(signal,1)/M), size(signal,2));
    for ch = 1:size(signal,2)
        y_decimated(:,ch) = decimate(signal(:,ch), M);
    end
else
    y_decimated = decimate(signal, M);
end
audiowrite('modulator22_downsampled_decimate.wav', y_decimated, round(Fs_ds));

% Temporal comparison plots (first channel if multi)
figure;
subplot(3,1,1);
plot((0:length(signal)-1)/Fs, signal(:,1));
title('Original (channel 1)');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,2);
plot((0:length(y_downsampled)-1)/Fs_ds, y_downsampled(:,1));
title('Downsample (downsample())');
xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,3);
plot((0:length(y_decimated)-1)/Fs_ds, y_decimated(:,1));
title('Decimate (decimate())');
xlabel('Time (s)'); ylabel('Amplitude');

% Q: Differences & distortions (store as comments for report)
% - downsample() simply keeps every M-th sample (no filtering) -> causes aliasing (spectral folding), possible high-frequency distortion.
% - decimate() applies an anti-aliasing IIR filter before downsampling -> reduces aliasing, smoother result but may introduce phase distortion.

% 3) Design low-pass FIR filter (order 30) with cutoff 1000 Hz
fir_order = 30;
fcut = 1000; % Hz
Wn = fcut / (Fs/2); % normalized cutoff (0..1)
b_fir = fir1(fir_order, Wn, 'low', hamming(fir_order+1));
a_fir = 1;

% 4) Design low-pass IIR Butterworth filter (order 8) with cutoff 1000 Hz
iir_order = 8;
[b_iir, a_iir] = butter(iir_order, Wn, 'low');

% Check stability:
% For FIR: all poles at origin (a_fir = 1) -> stable
% For IIR: check poles magnitude
p_iir = roots(a_iir);
stable_iir = all(abs(p_iir) < 1); % boolean

% Display stability results in command window
fprintf('FIR filter: order %d, stable (all poles at origin).\n', fir_order);
fprintf('IIR Butterworth: order %d, stable = %d (max pole magnitude = %.6f).\n', iir_order, stable_iir, max(abs(p_iir)));

% Display frequency responses (magnitude and phase)
figure;
freqz(b_fir, a_fir, 1024, Fs);
title('FIR Filter Frequency Response (mag & phase)');
figure;
freqz(b_iir, a_iir, 1024, Fs);
title('IIR Butterworth Frequency Response (mag & phase)');

% Apply FIR and IIR filters to original signal (filter each channel)
if size(signal,2) > 1
    y_fir = zeros(size(signal));
    y_iir = zeros(size(signal));
    for ch = 1:size(signal,2)
        y_fir(:,ch) = filter(b_fir, a_fir, signal(:,ch));
        y_iir(:,ch) = filter(b_iir, a_iir, signal(:,ch));
    end
else
    y_fir = filter(b_fir, a_fir, signal);
    y_iir = filter(b_iir, a_iir, signal);
end

audiowrite('modulator22_fir_filtered.wav', y_fir, Fs);
audiowrite('modulator22_iir_filtered.wav', y_iir, Fs);

% Compare temporal variation (first channel)
figure;
subplot(3,1,1); plot((0:length(signal)-1)/Fs, signal(:,1)); title('Original (ch1)');
subplot(3,1,2); plot((0:length(y_fir)-1)/Fs, y_fir(:,1)); title('FIR filtered (ch1)');
subplot(3,1,3); plot((0:length(y_iir)-1)/Fs, y_iir(:,1)); title('IIR filtered (ch1)');

% 5) Downsample the filtered signals to 4000 Hz using downsample()
y_fir_ds = downsample(y_fir, M);
y_iir_ds = downsample(y_iir, M);
audiowrite('modulator22_fir_downsampled.wav', y_fir_ds, round(Fs_ds));
audiowrite('modulator22_iir_downsampled.wav', y_iir_ds, round(Fs_ds));

% Plot temporal comparisons after filtering+downsampling (first channel)
figure;
subplot(3,1,1);
plot((0:length(y_downsampled)-1)/Fs_ds, y_downsampled(:,1)); title('Direct downsample() result');
subplot(3,1,2);
plot((0:length(y_fir_ds)-1)/Fs_ds, y_fir_ds(:,1)); title('FIR -> downsample()');
subplot(3,1,3);
plot((0:length(y_iir_ds)-1)/Fs_ds, y_iir_ds(:,1)); title('IIR -> downsample()');

