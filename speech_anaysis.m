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
