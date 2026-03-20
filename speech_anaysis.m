
% =========================================================
% 4.2 - Speech Analysis
% =========================================================

[signal, Fs] = audioread('modulator22.wav');

% If stereo, take first channel only
if size(signal, 2) > 1
    signal = signal(:, 1);
end

% ---------------------------------------------------------
% Save with half sampling frequency (same samples, new Fs)
% Effect: signal plays at half speed / half pitch
% ---------------------------------------------------------
Fs_half = round(Fs / 2);
audiowrite('modulator22_halfFs.wav', signal, Fs_half);

% Q: Halving Fs while keeping samples unchanged makes the player interpret
%    the signal as slower -> audio plays at half speed and half pitch.

% ---------------------------------------------------------
% 4.2.1 Temporal analysis
% ---------------------------------------------------------
t = (0:length(signal)-1).' / Fs;

figure;
plot(t, signal);
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');
title('Temporal variations of speech signal amplitude');
grid on;

% Q: The signal is NOT stationary. Spectral characteristics change over time
%    (vowels vs consonants vs silence). Speech is quasi-stationary over ~20-30 ms.

% ---------------------------------------------------------
% 4.2.2 Spectral analysis - DFT (pad=false)
% ---------------------------------------------------------
[power_dft, duration_dft] = frequencySpectrum(signal, Fs, false);
sgtitle('Power Spectrum - DFT (no padding)');

% Q: From the power spectrum we can retrieve: dominant frequencies,
%    fundamental frequency F0, formant peaks, and overall spectral shape.

% ---------------------------------------------------------
% Spectral analysis - FFT (pad=true, zero-pad to next power of 2)
% ---------------------------------------------------------
[power_fft, duration_fft] = frequencySpectrum(signal, Fs, true);
sgtitle('Power Spectrum - FFT (zero-padded to next power of 2)');

% Q: Setting pad=true zero-pads the signal to the next power of 2, enabling
%    the Cooley-Tukey FFT algorithm: O(N log N) instead of O(N^2) for DFT.

% Repeat timing measurements over 5 trials for statistics
ntrials = 5;
n_dft   = length(signal);
n_fft   = 2^nextpow2(n_dft);

dftTimes = zeros(ntrials, 1);
fftTimes = zeros(ntrials, 1);

for k = 1:ntrials
    tic; fft(signal, n_dft); dftTimes(k) = toc;
    tic; fft(signal, n_fft); fftTimes(k) = toc;
end

dft_mean = mean(dftTimes); dft_std = std(dftTimes);
fft_mean = mean(fftTimes); fft_std = std(fftTimes);

fprintf('DFT: mean = %.6f s, std = %.6f s\n', dft_mean, dft_std);
fprintf('FFT: mean = %.6f s, std = %.6f s\n', fft_mean, fft_std);

% Q: Yes, FFT is faster. DFT complexity is O(N^2); FFT (Cooley-Tukey) is
%    O(N log N), which is substantially faster for large N.

% ---------------------------------------------------------
% Spectrogram - window = 30 ms, step = 5 ms
% ---------------------------------------------------------
win_ms   = 30;  step_ms = 5;
win_samp = round(win_ms  / 1000 * Fs);
noverlap = win_samp - round(step_ms / 1000 * Fs);
nfft_s   = max(256, 2^nextpow2(win_samp));

figure;
spectrogram(signal, hamming(win_samp), noverlap, nfft_s, Fs, 'yaxis');
title(sprintf('Spectrogram (window=%d ms, step=%d ms)', win_ms, step_ms));
colorbar;

% Q: A spectrogram is a time-frequency-energy representation computed via
%    STFT: the signal is split into short overlapping windows, the FFT of
%    each window gives spectral content at that time instant. The result
%    shows how the spectrum evolves over time.
%    Difference from DFT: DFT gives one global spectrum; spectrogram tracks
%    how spectral content changes over time.

% ---------------------------------------------------------
% Spectrogram - window = 5 ms, step = 5 ms
% ---------------------------------------------------------
win_ms2   = 5;
win_samp2 = round(win_ms2 / 1000 * Fs);
noverlap2 = 0;  % step = window size, no overlap
nfft_s2   = max(256, 2^nextpow2(win_samp2));

figure;
spectrogram(signal, hamming(win_samp2), noverlap2, nfft_s2, Fs, 'yaxis');
title(sprintf('Spectrogram (window=%d ms, step=%d ms)', win_ms2, step_ms));
colorbar;

% Q: Shorter window -> better time resolution but worse frequency resolution
%    (uncertainty principle: Δt·Δf ≥ constant).
%    With 5 ms window, frequency bins are wide, formants blur together.
%    Optimal for speech: ~20-30 ms windows balance time and freq resolution.

% ---------------------------------------------------------
% Vowel analysis: /ʌ/ (one), /uː/ (two), /iː/ (three)
% !!! Adjust vowel_ranges (in seconds) after inspecting temporal plot !!!
% ---------------------------------------------------------
vowel_names  = {'/\Lambda/ (one)', '/u:/ (two)', '/i:/ (three)'};
vowel_ranges = [   % [start_s, end_s] - tune to your file
    0.25, 0.40;    % /ʌ/ in "one"
    0.65, 0.85;    % /uː/ in "two"
    1.05, 1.20;    % /iː/ in "three"
];

figure;
for v = 1:3
    s1  = max(1, round(vowel_ranges(v,1) * Fs) + 1);
    s2  = min(length(signal), round(vowel_ranges(v,2) * Fs));
    seg = signal(s1:s2);
    N   = length(seg);
    Np  = 2^nextpow2(N);

    Y  = fft(seg, Np);
    P2 = (abs(Y) / N).^2;
    P1 = P2(1:floor(Np/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    fv = (0:floor(Np/2)).' * (Fs / Np);

    % Smooth for formant envelope detection (~50 Hz smoothing)
    sw = max(3, round(50 / (Fs / Np)));
    P1s = movmean(P1, sw);

    % F0: dominant peak 50-500 Hz
    m0 = fv > 50 & fv < 500;
    [~, r] = max(P1(m0));
    F0 = fv(m0); F0 = F0(r);

    % F1, F2: top 2 local maxima of smooth spectrum 200-4000 Hz
    mf = fv > 200 & fv < 4000;
    Pf = P1s(mf); ff = fv(mf);
    pk = [false; Pf(2:end-1) > Pf(1:end-2) & Pf(2:end-1) > Pf(3:end); false];
    pf = ff(pk); pp = Pf(pk);
    [~, ord] = sort(pp, 'descend');
    tops = sort(pf(ord(1:min(2,numel(ord)))));
    F1 = tops(1);
    F2 = NaN; if numel(tops) >= 2, F2 = tops(2); end

    fprintf('Vowel %s: F0=%.1f Hz  F1=%.1f Hz  F2=%.1f Hz\n', vowel_names{v}, F0, F1, F2);

    subplot(1, 3, v);
    plot(fv, 10*log10(P1+eps), 'b', 'LineWidth', 0.5); hold on;
    plot(fv, 10*log10(P1s+eps), 'r', 'LineWidth', 1.5);
    xline(F0, 'g--', sprintf('F0=%.0f',F0));
    xline(F1, 'm--', sprintf('F1=%.0f',F1));
    if ~isnan(F2), xline(F2, 'k--', sprintf('F2=%.0f',F2)); end
    xlim([0 4000]);
    xlabel('Frequency (Hz)'); ylabel('Power (dB)');
    title(sprintf('Spectrum %s', vowel_names{v}));
    legend('Raw','Envelope','Location','NE'); grid on;
end
sgtitle('Vowel power spectra: F0, F1, F2 estimation');

% Reference averages (English):  /ʌ/ F1~700 F2~1200 Hz
%                                 /uː/ F1~300 F2~800 Hz
%                                 /iː/ F1~300 F2~2300 Hz

% ---------------------------------------------------------
% 4.2.3 Downsampling to 4000 Hz
% ---------------------------------------------------------
Fs_target = 4000;
M         = round(Fs / Fs_target);
Fs_ds     = Fs / M;

% Method 1: downsample() - no anti-aliasing filter
y_ds  = downsample(signal, M);
audiowrite('modulator22_downsampled_downsample.wav', y_ds, round(Fs_ds));

% Method 2: decimate() - low-pass filters before decimation
y_dec = decimate(signal, M);
audiowrite('modulator22_downsampled_decimate.wav', y_dec, round(Fs_ds));

figure;
subplot(3,1,1); plot((0:length(signal)-1)/Fs,   signal); title('Original');     xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,2); plot((0:length(y_ds)-1)/Fs_ds,  y_ds);   title('downsample()'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,3); plot((0:length(y_dec)-1)/Fs_ds, y_dec);  title('decimate()');   xlabel('Time (s)'); ylabel('Amplitude');
sgtitle('Downsampling: downsample() vs decimate()');

% Q: downsample() keeps every M-th sample with NO filtering -> aliasing.
%    decimate() applies an anti-aliasing low-pass filter first -> clean result.
%    Distortions from downsample(): high-frequency components above Nyquist
%    (2000 Hz) fold back into the spectrum (aliasing / spectral folding).

% ---------------------------------------------------------
% FIR filter: order 30, fc = 1000 Hz
% ---------------------------------------------------------
fir_order = 30;
f_cut     = 1000;
Wn        = f_cut / (Fs / 2);   % normalized [0..1]

b_fir = fir1(fir_order, Wn, 'low', hamming(fir_order + 1));
a_fir = 1;

% IIR Butterworth: order 8, fc = 1000 Hz
iir_order = 8;
[b_iir, a_iir] = butter(iir_order, Wn, 'low');

% Stability check
% FIR: denominator = 1 -> all poles at z=0 -> always stable
p_iir       = roots(a_iir);
stable_iir  = all(abs(p_iir) < 1);
fprintf('FIR  (order %d): always stable (all poles at origin).\n', fir_order);
fprintf('IIR  (order %d): stable = %d, max|pole| = %.6f\n', iir_order, stable_iir, max(abs(p_iir)));

% Q: Stability test: all poles must lie strictly inside the unit circle
%    (|pole| < 1) in the z-plane. FIR filters are unconditionally stable.
%    This Butterworth IIR should also be stable by design.

% Frequency responses
figure; freqz(b_fir, a_fir, 1024, Fs);
title(sprintf('FIR frequency response (order=%d, fc=%d Hz)', fir_order, f_cut));
figure; freqz(b_iir, a_iir, 1024, Fs);
title(sprintf('IIR Butterworth frequency response (order=%d, fc=%d Hz)', iir_order, f_cut));

% Apply filters
y_fir = filter(b_fir, a_fir, signal);
y_iir = filter(b_iir, a_iir, signal);
audiowrite('modulator22_fir_filtered.wav', y_fir, Fs);
audiowrite('modulator22_iir_filtered.wav', y_iir, Fs);

figure;
subplot(3,1,1); plot(t, signal); title('Original');     xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,2); plot(t, y_fir);  title('FIR filtered'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,3); plot(t, y_iir);  title('IIR filtered'); xlabel('Time (s)'); ylabel('Amplitude');
sgtitle('FIR vs IIR filtered signals (fc=1000 Hz)');

% Downsample filtered signals to 4000 Hz
y_fir_ds = downsample(y_fir, M);
y_iir_ds = downsample(y_iir, M);
audiowrite('modulator22_fir_downsampled.wav', y_fir_ds, round(Fs_ds));
audiowrite('modulator22_iir_downsampled.wav', y_iir_ds, round(Fs_ds));

figure;
subplot(3,1,1); plot((0:length(y_ds)-1)/Fs_ds,     y_ds);     title('Direct downsample() only');  xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,2); plot((0:length(y_fir_ds)-1)/Fs_ds, y_fir_ds); title('FIR filter + downsample()'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,3); plot((0:length(y_iir_ds)-1)/Fs_ds, y_iir_ds); title('IIR filter + downsample()'); xlabel('Time (s)'); ylabel('Amplitude');
sgtitle('Effect of anti-aliasing filter before downsampling');

% Q: Conclusion: pre-filtering with a low-pass filter (FIR or IIR) before
%    downsampling eliminates aliasing. IIR achieves sharper roll-off at lower
%    order but introduces non-linear phase. FIR has linear phase but needs
%    higher order. Both produce cleaner downsampled audio than raw downsample().