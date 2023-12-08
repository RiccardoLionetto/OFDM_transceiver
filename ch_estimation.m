function h_hat = ch_estimation(rxsignal, conf)

% Incoming signal is composed as follows: [CP TRAINING CP DATA]

training_seq = lfsr_training(conf.train_length);
training_symb = -1 +2*training_seq;

time = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
rx_downmixing = rxsignal.*exp(-2*1j*pi*conf.f_c*time);

% Remove CP []
rxxx = rx_downmixing(conf.Ncp+1:end);
%filtered_signal = ofdmlowpass(rxxx,conf,conf.carriers*conf.spacing); % 256*5
%rx_symbols = osfft(filtered_signal, conf.os_factor);
rx_symbols = osfft(rxxx, conf.os_factor);

rx_symbols = rx_symbols*2;

% rx_noCP= CP_OFDM_symb(conf.Ncp+1:end);
% rx_symbols = osfft(rx_noCP, conf.os_factor);


figure
plot(rx_symbols, '.');
title("rx symbols after osfft")

h_hat = rx_symbols(1:conf.train_length)./training_symb;




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
