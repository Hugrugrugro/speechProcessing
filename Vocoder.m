
% =========================================================
% 4.3 - Channel Vocoder
% Download chanvocoder.m from:
% https://sethares.engr.wisc.edu/vocoders/channelvocoder.html
% Place chanvocoder.m in the same directory as this script.
% =========================================================

% ---------------------------------------------------------
% Parameters (adjust to experiment)
% ---------------------------------------------------------
num_channels = 64;   % half the FFT length -> frequency resolution
num_bands    = 16;   % number of bandpass filter bands
win_overlap  = 0.5;  % window overlap fraction (0..1)

% ---------------------------------------------------------
% Test 1: modulator = voice, carrier = synthesizer
% ---------------------------------------------------------
modulator_file = 'modulator22.wav';
carrier_file   = 'carrier22.wav';  % provided synthesizer carrier

[modulator, Fs_mod] = audioread(modulator_file);
[carrier,   Fs_car] = audioread(carrier_file);

% Ensure mono
if size(modulator, 2) > 1, modulator = modulator(:,1); end
if size(carrier,   2) > 1, carrier   = carrier(:,1);   end

% Resample carrier to match modulator Fs if needed
if Fs_car ~= Fs_mod
    carrier = resample(carrier, Fs_mod, Fs_car);
end

% Match lengths (truncate to shorter)
min_len   = min(length(modulator), length(carrier));
modulator = modulator(1:min_len);
carrier   = carrier(1:min_len);

% Apply channel vocoder
output_synth = chanvocoder(modulator, carrier, num_channels, num_bands, win_overlap, Fs_mod);
audiowrite('vocoder_output_synth.wav', output_synth / max(abs(output_synth)), Fs_mod);

% Temporal comparison plot
t = (0:min_len-1).' / Fs_mod;
figure;
subplot(3,1,1); plot(t, modulator);    title('Modulator (voice)');      xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,2); plot(t, carrier);      title('Carrier (synthesizer)');  xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,3); plot(t, output_synth); title('Vocoder output');         xlabel('Time (s)'); ylabel('Amplitude');
sgtitle('Channel Vocoder - synthesizer carrier');

% Spectrogram comparison
win_ms      = 30;
win_samp    = round(win_ms / 1000 * Fs_mod);
noverlap    = round(win_samp * 0.75);
nfft_spec   = max(256, 2^nextpow2(win_samp));

figure;
subplot(1,3,1); spectrogram(modulator,    hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Modulator');       colorbar;
subplot(1,3,2); spectrogram(carrier,      hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Carrier');         colorbar;
subplot(1,3,3); spectrogram(output_synth, hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Vocoder output');  colorbar;
sgtitle('Spectrograms: modulator / carrier / vocoder output (synthesizer)');

% ---------------------------------------------------------
% Test 2: white noise carrier
% ---------------------------------------------------------
white_noise = randn(min_len, 1);
white_noise = white_noise / max(abs(white_noise));

output_noise = chanvocoder(modulator, white_noise, num_channels, num_bands, win_overlap, Fs_mod);
audiowrite('vocoder_output_noise.wav', output_noise / max(abs(output_noise)), Fs_mod);

figure;
subplot(1,3,1); spectrogram(modulator,    hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Modulator');           colorbar;
subplot(1,3,2); spectrogram(white_noise,  hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('White noise carrier'); colorbar;
subplot(1,3,3); spectrogram(output_noise, hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Vocoder output');      colorbar;
sgtitle('Spectrograms: modulator / carrier / vocoder output (white noise)');

% ---------------------------------------------------------
% Test 3: periodic (pulse train) noise carrier -> more "robotic"
% ---------------------------------------------------------
period_samp    = round(Fs_mod / 120);  % 120 Hz fundamental (robotic pitch)
periodic_noise = zeros(min_len, 1);
periodic_noise(1:period_samp:end) = 1;
% Smooth with narrow Gaussian to reduce click harshness
periodic_noise = conv(periodic_noise, gausswin(15)/sum(gausswin(15)), 'same');
periodic_noise = periodic_noise / max(abs(periodic_noise));

output_periodic = chanvocoder(modulator, periodic_noise, num_channels, num_bands, win_overlap, Fs_mod);
audiowrite('vocoder_output_periodic.wav', output_periodic / max(abs(output_periodic)), Fs_mod);

figure;
subplot(1,3,1); spectrogram(modulator,       hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Modulator');              colorbar;
subplot(1,3,2); spectrogram(periodic_noise,  hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Periodic noise carrier'); colorbar;
subplot(1,3,3); spectrogram(output_periodic, hamming(win_samp), noverlap, nfft_spec, Fs_mod, 'yaxis'); title('Vocoder output');         colorbar;
sgtitle('Spectrograms: modulator / carrier / vocoder output (periodic noise)');

% Q: The vocoder extracts the spectral envelope (formant structure) from the
%    modulator and imposes it on the carrier's energy. The output retains the
%    speech intelligibility of the modulator with the timbre/texture of the carrier.
%    - Synthesizer carrier: robotic but tonal, pitched output.
%    - White noise carrier:  whispery/noisy speech (all voiced sounds become unvoiced).
%    - Periodic noise carrier: robotic monotone, fixed F0 regardless of original pitch.

% Q: Effect of parameters:
%    - num_channels: more channels = finer frequency resolution = more spectral detail.
%    - num_bands: more bands = more granular envelope shaping across frequency.
%    - win_overlap: higher overlap = smoother temporal transitions in the output.

% ---------------------------------------------------------
% Live microphone input -> vocoder -> playback
% (Requires audio hardware; comment out if not available)
% ---------------------------------------------------------
rec_duration = 3;   % seconds to record
fprintf('Recording %d seconds from microphone...\n', rec_duration);

try
    rec    = audiorecorder(Fs_mod, 16, 1);  % Fs, bits, channels
    recordblocking(rec, rec_duration);
    mic_signal = getaudiodata(rec, 'double');

    % Trim/pad carrier to match recorded length
    mic_len = length(mic_signal);
    if length(carrier) >= mic_len
        carrier_live = carrier(1:mic_len);
    else
        carrier_live = [carrier; zeros(mic_len - length(carrier), 1)];
    end

    output_live = chanvocoder(mic_signal, carrier_live, num_channels, num_bands, win_overlap, Fs_mod);
    output_live = output_live / max(abs(output_live) + eps);

    fprintf('Playing vocoded output...\n');
    player = audioplayer(output_live, Fs_mod);
    playblocking(player);

    audiowrite('vocoder_live_output.wav', output_live, Fs_mod);
    fprintf('Saved to vocoder_live_output.wav\n');
catch e
    fprintf('Microphone recording failed: %s\n', e.message);
    fprintf('Skipping live recording section.\n');
end
