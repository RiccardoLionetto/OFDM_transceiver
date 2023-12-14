function txpreamble = tx_preamble(conf)

% Generatate preamble and map it to BPSK
preamble_mapped = -1 + 2*lfsr_framesync(conf.npreamble);

% Up sampling
preamble_upsampled = upsample(preamble_mapped,conf.os_factor_preamble);

% Pulse shaping filtering
pulse = rrc(conf.os_factor_preamble, conf.roll_off, conf.filterlength);
filtered_txsignal = conv(preamble_upsampled, pulse, 'same');

% figure
% plot(abs(fftshift(fft(filtered_txsignal))));
% title("pulse-shaped preamble, in freq domain")
% 
% figure
% plot(filtered_txsignal)
% title("pulse-shaped preamble, in time domain")

%UP CONVERSION
time = (0:1/conf.f_s:(length(filtered_txsignal)-1)/conf.f_s)';
txpreamble = real(filtered_txsignal).*cos(2*pi*conf.f_c.*time)-imag(filtered_txsignal).*sin(2*pi*conf.f_c.*time);

% Normalizing signal with respect to average energy
% average_energy = norm(txpreamble)^2 /length(txpreamble);
% txpreamble = txpreamble/sqrt(average_energy);

% figure
% plot(abs(fftshift(fft(txpreamble))));
% title("up-converted preamble")

end