function training_txsignal = training_frame(conf)

% Generation of sequence
training_seq = lfsr_training(conf.train_length);

% Mapping in BPSK
symbols = -1 +2*training_seq;

% Over-sampled IFFT
OFDM_symbols = osifft(symbols,conf.os_factor);

% Add Cyclic Prefix
CP_OFDM_symb = zeros(1, length(OFDM_symbols)+conf.Ncp);
CP = OFDM_symbols(end-conf.Ncp+1: end);
CP_OFDM_symb(1:conf.Ncp) = CP;
CP_OFDM_symb(conf.Ncp+1:end) = OFDM_symbols;

% Mixing
time = (0:1/conf.f_s:(length(CP_OFDM_symb)-1)/conf.f_s);
%training_txsignal = real(CP_OFDM_symb(1:end)).*cos(2*pi*conf.f_c.*time(1:length(CP_OFDM_symb)))-imag(CP_OFDM_symb(1:end)).*sin(2*pi*conf.f_c.*time(1:length(CP_OFDM_symb)));
training_txsignal = real(CP_OFDM_symb(1:end)).*cos(2*pi*conf.f_c.*time)-imag(CP_OFDM_symb(1:end)).*sin(2*pi*conf.f_c.*time);



%to cancel
%rxsignal = training_txsignal;

%time = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
% RX

%rx_downmixing = rxsignal.*exp(-2*1j*pi*conf.f_c*time);
% Remove CP []
%
%rxxx = rx_downmixing(conf.Ncp+1:end);
%filtered_signal = ofdmlowpass(rxxx,conf,conf.carriers*conf.spacing); % 256*5
%rx_symbols = osfft(filtered_signal, conf.os_factor);
%rx_symbols = osfft(rxxx, conf.os_factor);

%rx_symbols = rx_symbols*2;

% rx_noCP= CP_OFDM_symb(conf.Ncp+1:end);
% rx_symbols = osfft(rx_noCP, conf.os_factor);


%figure
%plot(rx_symbols, '.');
%title("rx symbols after osfft")

%h_hat = rx_symbols(1:conf.train_length)./symbols;



end