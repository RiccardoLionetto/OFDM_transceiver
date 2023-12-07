function [rxbits conf] = rx_ofdm(rxsignal,conf, h)

% remove 
%rxsignal = rxsignal(length(rxsignal)/2+1:end);
% Down mixing
time_rx = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
rx_downmixing = rxsignal.*exp(-2*1j*pi*conf.f_c*time_rx).';
% Remove CP
rxsignal = rx_downmixing(conf.Ncp+1:end);
filtered_signal = ofdmlowpass(rxsignal,conf,conf.carriers*conf.spacing); % 256*5

% Transform
rx_symbols = osfft(filtered_signal, conf.os_factor);

% Channel correction: DOESN'T WORK
% rx_symbols = rx_symbols./h;

figure
plot(rx_symbols, '.');
title("rx symbols after osfft")

% Demapping QPSK
constellation = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
[~,ind] = min(abs(ones(conf.bitsXsymb/2,4)*diag(constellation) - diag(rx_symbols)*ones(conf.bitsXsymb/2,4)),[],2);
rxbits = de2bi(ind-1, 'left-msb',2);
rxbits = rxbits(:);

% BER
biterrors = sum(rxbits ~= txbits);




% 
% rxsignal_noCP = rxsignal(conf.Ncp:end);
% 
% % Matched filtering
% rxsignal_baseband = rxsignal_noCP.*exp(-2*1j*pi*conf.f_c*(0:1/conf.f_s:(length(rxsignal_noCP)-1)/conf.f_s).');
% rxsignal = ofdmlowpass(rxsignal_noCP.*exp(-2*1j*pi*conf.f_c*(0:1/conf.f_s:(length(rxsignal_noCP)-1)/conf.f_s).'), conf, 256);
% 
% figure
% plot(abs(fft(rxsignal_baseband)));
% title("freq domain rx signal baseband")
% 
% figure
% plot(abs(fft(rxsignal)));
% title("freq domain rx signal baseband lowpass")


end