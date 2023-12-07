function [txsignal conf] = tx_ofdm(txbits,conf)

if conf.mapping == 4 % this will be called when generating data
    % Map bits in QPSK
    symbols = QPSK(txbits); % Vector of dimensions 256x1 in complex domain
elseif conf.mapping == 2 % this will be called when generating training
    % Map bits in BPSK
    symbols = BPSK(txbits);
else
    disp('ERROR: mapping is possible only on BPSK and QPSK'); 
end

% Over-sampled IFFT
OFDM_symbols = osifft(symbols,conf.os_factor);

% Add Cyclic Prefix
CP_OFDM_symb = zeros(1, length(OFDM_symbols)+conf.Ncp);
CP = OFDM_symbols(end-conf.Ncp+1: end);
CP_OFDM_symb(1:conf.Ncp) = CP;
CP_OFDM_symb(conf.Ncp+1:end) = OFDM_symbols;

% Mixing
time = (0:1/conf.f_s:(length(CP_OFDM_symb)-1)/conf.f_s);
conf.time = time;
txsignal = real(CP_OFDM_symb).*cos(2*pi*conf.f_c.*time)-imag(CP_OFDM_symb).*sin(2*pi*conf.f_c.*time);

figure
plot(time, txsignal)
title("sent signal in time")

% figure
% plot(1:5:9600*5, abs(fft(txsignal)));
% title("tx in freq")


% % RX
% rxsignal = txsignal;
% % Down mixing
% rx_downmixing = rxsignal.*exp(-2*1j*pi*conf.f_c*conf.time);
% % Remove CP
% filtered_signal = ofdmlowpass(rx_downmixing,conf,conf.carriers*conf.spacing*conf.os_factor); % 256*5
% rxsignal = filtered_signal(conf.Ncp+1:end);
% 
% 
% % Transform
% rx_symbols = osfft(rxsignal, conf.os_factor);
% 
% figure
% plot(rx_symbols, '.');
% title("rx symbols after osfft")
% 
% % Demapping QPSK
% constellation = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
% [~,ind] = min(abs(ones(conf.bitsXsymb/2,4)*diag(constellation) - diag(rx_symbols)*ones(conf.bitsXsymb/2,4)),[],2);
% rxbits = de2bi(ind-1, 'left-msb',2);
% rxbits = rxbits(:);
% 
% % BER
% biterrors = sum(rxbits ~= txbits);

% ------------------------------------------------------------
% figure
% plot(symbols,'x');
% title("mapped symbols - expecting QPSK")
% 
% 
% % IFFT
% %OFDM_symbols = osifft(symbols,conf.os_factor);
% OFDM_symbols = ifft(symbols);
% 
% % Definition of txsignal: it will be composed by Cyclic Prefix followed by
% % the OFDM symbol
% txsignal = zeros(1, length(OFDM_symbols)+conf.Ncp);
% 
% figure
% plot(OFDM_symbols,'.')
% title("time domain after IFFT")
% 
% figure
% plot(abs(OFDM_symbols))
% title("time domain after IFFT")
% 
% % Add Cyclic Prefix
% %CP = OFDM_symbols(end-conf.Ncp*conf.f_s+2: end);
% CP = OFDM_symbols(end-conf.Ncp+1: end);
% %txsignal(1:conf.Ncp) = CP;
% txsignal(conf.Ncp+1:end) = OFDM_symbols;
% 
% % figure
% % plot(txsignal,'.')
% % title("time domain after add of CP")
% 
% 
% % Mixing
% time = (0:length(txsignal)-1);
% txsignal = real(txsignal(1:end)).*cos(2*pi*conf.f_c.*time(1:length(txsignal)))-imag(txsignal(1:end)).*sin(2*pi*conf.f_c.*time(1:length(txsignal)));
% 
% 
% figure
% plot(time, txsignal)
% title("time domain tx signal final")
% 
% figure
% plot(abs(fft(txsignal)));
% title("freq domain tx signal final")
% 
% 
% % receving - to be deleted
% rx_downmixing = txsignal.*exp(-2*1j*pi*conf.f_c*time);
% filtered_signal = ofdmlowpass(rx_downmixing,conf,400);
% 
% 
% figure
% plot(abs(fft(rx_downmixing)));
% title("rx_downmixing")
% 
% figure
% plot(abs(fft(filtered_signal)));
% title("filtered_signal")
% 
% figure
% plot(fft(filtered_signal));
% title("received symbols filtered",'x')
% 
% 
% figure
% plot(fft(rx_downmixing));
% title("received symbols not filtered",'x')

% --------------------------------------- %
% Definition of functions

% QPSK function mapping
function[symbol] = QPSK (txbits)
    txbits = reshape(txbits,[],2);
    GrayMap = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
    symbol = GrayMap(bi2de(txbits, 'left-msb')+1).';
end

function[symbol] = BPSK (txbits)
    symbol = -1 + 2*txbits;
end



end