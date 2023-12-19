function beginning_of_data = frame_sync(rxsignal, conf)

% Frame synchronizer: it is used to only detect the beginning of the
% sequence. .
%% Useful plots
% Processing on rx: down-mixing, lowplass, matched filtering
% time = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s)';
% 
% figure
% plot(abs(fftshift(fft(rxsignal))));
% title("rx signal  - freq domain")
% figure
% plot(rxsignal);
% title("rx signal - time domain")

% figure
% plot(abs(fftshift(fft(rxsignal))));
% title("rx signal after downmixing  - freq domain ")
% figure
% plot(rxsignal);
% title("rx signal after downmixing - time domain")

% 
% figure
% plot(abs(fftshift(fft(rx_signal))));
% title("rx signal after pulse shaping - freq domain")
% 
% figure
% plot(rx_signal);
% title("rx signal after pulse shaping - time domain")

%% 
% Exclude the preamble from the carrier frequency offset
time_rx = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
rxsignal = rxsignal .*exp(+2*1j*pi*(conf.carrierOffset)*time_rx);

rxsignal = lowpass_preamble(rxsignal.', conf);
pulse = rrc(conf.os_factor_preamble, conf.roll_off, conf.filterlength);
rx_signal = conv(rxsignal, pulse, 'same');

% Define values for detection
detection_threshold = 50;
L = conf.os_factor_preamble;
% Regenerate same preamble of tx
frame_sync_sequence = -1 + 2*lfsr_framesync(conf.npreamble);
frame_sync_length = conf.npreamble;
% Initialize support variables
current_peak_value = 0;
samples_after_threshold = L;
vector_T = zeros(1,frame_sync_length);


for i = L * frame_sync_length + 1 : length(rx_signal)
    r = rx_signal(i - L * frame_sync_length : L : i - L);
    c = frame_sync_sequence' * r;
    T = abs(c)^2 / abs(r' * r);
    vector_T(i)=T;
    
    if (T > detection_threshold || samples_after_threshold < L)
        samples_after_threshold = samples_after_threshold - 1;
        if (T > current_peak_value)
            beginning_of_data = i;
            current_peak_value = T;
        end
        if (samples_after_threshold == 0)
            figure
            plot(1:length(vector_T), vector_T)
            title("vector T - frame synch - it worked")
            return;
        end
    end
end
figure
plot(1:length(vector_T), vector_T)
title("vector T - frame synch - didn't get it")

error('No synchronization sequence found.');
return