function [power, duration] = frequencySpectrum(signal, fs, pad)
%%%%%%%%%%%%%%%%%%
%function power = frequencySpectrum(signal, fs, pad)
%
% Task: Display the power spectrum (lin and log scale) of a given signal
%
% Input:
%	- signal: the input signal to process
%	- fs: the sampling rate
%	-pad: boolean if true, signal is padded with 0 to the next power of 2 -> FFT instead of DFT
%
% Output: 
%	- power: the power spectrum
%	
%
% Guillaume Gibert, guillaume.gibert@ecam.fr
% 25/04/2022
%%%%%%%%%%%%%%%%%%

n = length(signal);        % number of samples

if (pad)
	n = 2^nextpow2(n);
end

tic
y = fft(signal, n);% compute DFT of input signal
duration = toc;

power = abs(y).^2/n;    % power of the DFT

[val, ind] = max(power); % find the mx value of DFT and its index

% plots
figure;

subplot(1,3,1) % time plot 
t=0:1/fs:(n-1)/fs; % time range
%pad signal with zeros
if (pad)
	signal = [ signal; zeros( n-length(signal), 1)];
end
plot(t, signal)
xticks(0:0.1*fs:n*fs);
xticklabels(0:0.1:n/fs);
xlabel('Time (s)');
ylabel('Amplitude (a.u.)');

subplot(1,3,2) % linear frequency plot
f = (0:n-1)*(fs/n);     % frequency range
plot(f,power, 'b*'); hold on;
plot(f,power, 'r');
xlabel('Frequency (Hz)')
ylabel('Power (a.u.)')

subplot(1,3,3) % log frequency plot
plot(f,10*log10(power/power(ind)));
xlabel('Frequency (Hz)')
ylabel('Power (dB)')

