function y = chanvocoder(modulator, carrier, chan, numband, overlap, fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function y = chanvocoder(modulator, carrier, chan, numband, overlap, fs)
%
% Channel vocoder based on W. Sethares' algorithm.
% Imposes the spectral envelope of the modulator onto the carrier.
%
% Inputs:
%   modulator : modulator signal (voice), column vector
%   carrier   : carrier signal (synth/noise), column vector
%   chan      : number of FFT channels (FFT length = 2*chan)
%   numband   : number of frequency bands
%   overlap   : window overlap fraction (0 < overlap < 1), e.g. 0.5
%   fs        : sampling frequency (Hz)
%
% Output:
%   y         : vocoded output signal
%
% Based on: W.A. Sethares, channelvocoder.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- ensure column vectors ---
modulator = modulator(:);
carrier   = carrier(:);

% --- match lengths ---
N = min(length(modulator), length(carrier));
modulator = modulator(1:N);
carrier   = carrier(1:N);

fftlen   = 2 * chan;                        % FFT length
hopsize  = round(fftlen * (1 - overlap));   % hop size in samples
win      = hann(fftlen, 'periodic');        % analysis/synthesis window

% band edges: numband+1 linearly spaced bin indices from 1 to chan+1
band_edges = round(linspace(1, chan + 1, numband + 1));

% output accumulator
y       = zeros(N + fftlen, 1);
win_sum = zeros(N + fftlen, 1);   % for overlap-add normalisation

num_frames = floor((N - fftlen) / hopsize) + 1;

for frame = 1:num_frames
    idx = (frame - 1) * hopsize + 1;          % start sample of this frame
    idx_end = idx + fftlen - 1;

    % windowed frames
    mod_frame = modulator(idx:idx_end) .* win;
    car_frame = carrier(idx:idx_end)   .* win;

    % FFT
    MOD = fft(mod_frame, fftlen);
    CAR = fft(car_frame, fftlen);

    % magnitude spectra (one-sided: bins 1..chan+1)
    mod_mag = abs(MOD(1:chan+1));
    car_mag = abs(CAR(1:chan+1));

    % build gain spectrum: for each band, ratio of modulator to carrier energy
    gain = zeros(chan + 1, 1);
    for b = 1:numband
        b_start = band_edges(b);
        b_end   = band_edges(b + 1);
        bins    = b_start:b_end;

        mod_energy = sum(mod_mag(bins).^2) + eps;
        car_energy = sum(car_mag(bins).^2) + eps;

        band_gain         = sqrt(mod_energy / car_energy);
        gain(bins)        = band_gain;
    end

    % apply gain to carrier spectrum (full two-sided)
    gain_full          = [gain; flipud(gain(2:end-1))];
    OUT                = CAR .* gain_full;

    % inverse FFT -> time domain frame
    out_frame = real(ifft(OUT, fftlen)) .* win;

    % overlap-add
    y(idx:idx_end)       = y(idx:idx_end)       + out_frame;
    win_sum(idx:idx_end) = win_sum(idx:idx_end) + win.^2;
end

% normalise overlap-add (avoid division by zero)
win_sum(win_sum < 1e-8) = 1;
y = y ./ win_sum;

% trim to original length and normalise output
y = y(1:N);
peak = max(abs(y));
if peak > 0
    y = y / peak * 0.99;
end

end
