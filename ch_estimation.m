function h_hat = ch_estimation(rxsignal, conf)

% Incoming signal is composed as follows: [CP TRAINING CP DATA]

training_seq = lfsr_training(conf.train_length);
training_symb = -1 +2*training_seq;

time = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s)';

% Down mixing
rx_downmixing = rxsignal.*exp(-2*1j*pi*conf.f_c*time);

% Low pass
filtered_signal = ofdmlowpass(rx_downmixing,conf,conf.carriers*conf.spacing); % 256*5

% Remove CP
rx_symbols_noCP = filtered_signal(conf.Ncp+1:end);

% Over-sampled FFT
rx_symbols = osfft(rx_symbols_noCP, conf.os_factor);

rx_symbols = rx_symbols*2;


figure
plot(rx_symbols, '.');
title("rx symbols after osfft")

h_hat = rx_symbols(1:conf.train_length)./training_symb;

% Fare la fft di h_hat e ricavare il Power delay profile (da qui trarre le
% considerazioni su lunghezza CP)



%---------

% training_seq = lfsr_framesync(conf.train_length);
% training_symb = -1 +2*training_seq;
% 
% time = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
% % % RX
% 
% % Down mixing
% rx_downmixing = rxsignal.*exp(-2*1j*pi*conf.f_c*time);
% filtered_signal = ofdmlowpass(rxsignal,conf,conf.carriers*conf.spacing); % 256*5
% % Remove CP []
% rxsignal = rx_downmixing(conf.Ncp+1:end);
% 
% 
% rx_symbols = osfft(filtered_signal, conf.os_factor);
% 
% figure
% plot(rx_symbols, ".")
% 
% h_hat = rx_symbols(1:conf.train_length)./training_symb;


end
