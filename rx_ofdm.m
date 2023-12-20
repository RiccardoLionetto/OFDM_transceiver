function [rxbits conf] = rx_ofdm(rxsignal,conf)


% --------------- STAGE A ---------------
% Down-mixing and introduction of a carrier frequency offset (equal to the
% value of conf.carrierOffset). If set to zero, perfect matching between tx
% and rx oscillators frequency is assumed.
time_rx = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
rx_downmixing = rxsignal.*exp(-2*1j*pi*(conf.f_c+conf.carrierOffset)*time_rx);

figure
plot(abs(fftshift(fft(rx_downmixing))))
title("Down-mixed rx signal - freq domain")

% -------------------- STAGE B -----------------
% Frame synchronization to extract the starting index
idx_start = frame_sync(rx_downmixing,conf);
rx_OFDM = rx_downmixing(idx_start:end);

% -------------------- STAGE C -----------------
% Lowpass filtering on the OFDM sequence
filtered_signal = ofdmlowpass(rx_OFDM,conf,ceil((conf.carriers+1)/2)*conf.spacing);

figure
plot(abs(fftshift(fft(filtered_signal))))
title("Down-mixed and filtered rx signal - freq domain")

%Serial to parallel conversion
ofdm_frame= filtered_signal(1:conf.len_txsignal);
OFDM_td_matrix = reshape(ofdm_frame, conf.carriers*conf.os_factor+conf.Ncp, []).';


% -------------------- STAGE D -----------------
% Cyclic prefix removal
OFDM_td_noCP_matrix = OFDM_td_matrix(:,1+conf.Ncp:end);
OFDM_td_noCP_matrix = reshape(OFDM_td_noCP_matrix.', 1, []);
OFDM_td_noCP_matrix = OFDM_td_noCP_matrix.';

% -------------------- STAGE E -----------------
% Over-sampled Fast Fourier Transform to go back into frequency domain
mat_line = [];
for ii=0:length(OFDM_td_noCP_matrix)/(conf.carriers*conf.os_factor) -1
    mat = osfft(OFDM_td_noCP_matrix((ii*conf.carriers*conf.os_factor)+1:(ii*conf.carriers*conf.os_factor)+conf.carriers*conf.os_factor),conf.os_factor);%*80;
    mat_line = [mat_line; mat];
end
OFDM_fd_vector = mat_line;

figure
plot(OFDM_fd_vector, '.')
title("Rx signal after osfft [there are data, BPSK training symbols, zero pads]")

OFDM_fd_vector = reshape(OFDM_fd_vector, conf.carriers, []).';


% -------------------- STAGE F  -----------------
% CHANNEL ESTIMATION
% Extract training symbols
training_symbols = OFDM_fd_vector(1:conf.trainOccurrence+1:end, :); 
% Recreate original training symbols to make the comparison
sent_training_sequence = (-1 +2*lfsr_training(conf.train_length))';
sent_training_matrix = repmat(sent_training_sequence, [size(training_symbols,1), 1]);
% Division to obtain the channel estimation
H_hat = training_symbols ./ sent_training_matrix;

% Eliminate training symbols from OFDM_fd_vector to extract data from incoming T+D stream
OFDM_fd_vector(1:conf.trainOccurrence+1:end, :) = [];
data_symbols = OFDM_fd_vector;

figure
plot(data_symbols, '.')
title("Rx signal after osfft and removal of training")


% --------------- STAGE G: Channel correction (Training + Viterbi-Viterbi) ---------------
% Correct them by using channel estimation: every symbol has to be
% corrected by using the channel estimation of the last sent training
% sequence. Some operations are therefore performed first to allow this computation in a matrix form
expanded_ch_matrix = repelem(H_hat, conf.trainOccurrence, 1);
expanded_ch_matrix = expanded_ch_matrix(1:size(data_symbols,1),:);

% ------------------------------—------------------------------—
theta_hat = zeros(size(expanded_ch_matrix,1), conf.carriers);   % Estimate of the carrier phase
data_ch_corrected = zeros(size(data_symbols));
final_ch_estimation_matrix = zeros(size(theta_hat,1)+size(H_hat,1), conf.carriers);
ss=1; % index used to store ch estimations inside final_ch_estimation_matrix

if conf.ViterbiViterbi == 1
% Apply viterbi-viterbi algorithm
for train_index = 1:size(H_hat, 1) % Loop over received training symbols
    theta_hat_train = angle(H_hat(train_index,:));
    final_ch_estimation_matrix(ss,:) = abs(H_hat(train_index,:)).*exp(1j*theta_hat_train);
    ss = ss+1;
    for index=(train_index-1)*conf.trainOccurrence+1:(train_index-1)*conf.trainOccurrence+conf.trainOccurrence
        if index > size(expanded_ch_matrix, 1)
            continue
        end
        deltaTheta = (1/4*angle(-data_symbols(index,:).^4)).';% + pi/2*(-1:4);
        deltaTheta = repelem(deltaTheta, 1 ,6);
        rotation = repelem(pi/2*(-1:4), conf.carriers,1);
        deltaTheta_rotations = deltaTheta+rotation;

        % Pick values closer to previous estimation (no fast change approx)
        if index == (train_index-1)*conf.trainOccurrence+1
            prev_theta = theta_hat_train;
            prev_theta_rep = repelem(prev_theta.', 1,6);
        else
            prev_theta = theta_hat(index-1, :);
            prev_theta_rep = repelem(prev_theta.', 1,6);
        end

        % Find closest values and corresponding phase
        [~, ind] = min(abs(deltaTheta_rotations-prev_theta_rep), [], 2);
        idx = sub2ind(size(deltaTheta_rotations),(1:numel(ind)).',ind);
        theta = deltaTheta_rotations(idx);

        % Phase correction (and magnitude from prev. ch symbol)
        if index ~= (train_index-1)*conf.trainOccurrence+1
        theta_hat(index,:) = mod(0.3*theta.' + 0.7*prev_theta, 2*pi);
        data_ch_corrected(index,:) = data_symbols(index,:) .* exp(-1j * theta_hat(index,:)) ./abs(H_hat(train_index,:));
        else
            % The estimation is taken directly from the last training symbol
            theta_hat(index,:) = prev_theta;
            data_ch_corrected(index,:) = data_symbols(index,:) .* exp(-1j * theta_hat(index,:)) ./abs(H_hat(train_index,:));
        end
        final_ch_estimation_matrix(ss,:) = abs(H_hat(train_index,:)).*exp(1j*theta_hat(index,:));
        ss = ss+1;

    end
end


figure
hold on
plot(rad2deg(angle(data_symbols(:,1))))
plot(rad2deg(theta_hat(:,1)))
plot(rad2deg(angle(data_ch_corrected(:,1))))
plot(rad2deg(angle(expanded_ch_matrix(:,1))))
hold off
title("Phase of the 1° subcarrier throughout rx - Viterbi on")
xlabel("OFDM symbol")
ylabel("Degree")
legend('Data phase before correction', 'Estimated phase correction - Viterbi', 'Data phase after correction', 'Estimated phase correction - no Viterbi')

else
% Instead of applying Viterbi, only phase estimated from Training symbols is used
data_ch_corrected = data_symbols ./ expanded_ch_matrix;

figure
hold on
plot(angle(expanded_ch_matrix(:,1)))
plot(angle(data_symbols(:,1)))
hold off
title("Phase of the 1° subcarrier throughout rx - Viterbi off")
xlabel("OFDM symbol")
ylabel("Degree")
legend('Estimated phase correction - no Viterbi', 'Data phase before correction')

end
% ------------------------------—------------------------------—

data_stream = reshape(data_ch_corrected.', 1, []);
data_stream = data_stream(1:end-conf.npad);% remove pad added in transmission

figure
plot(data_stream, '.')
title("rx signal after:osfft, removal of training & zero pads, ch correction")

%Power delay profile of 1° subcarrier over different OFDM symbols recevied
h_hat_fft = fft(H_hat(1,:));
h_hat_power = abs(h_hat_fft).^2;
figure;
plot((0:1/conf.f_s:(length(h_hat_power)-1)/conf.f_s), h_hat_power.', 'Marker', 'none', 'LineWidth', 1.5);
title('Power Delay Profile (PDP)');
xlabel('Delay [second]');
ylabel('Power');

% -------------------- STAGE G ---------------
% Apply QPSK demapping over the data using Maximum-Likelihood estimation
constellation = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
demapping_matrix = repmat(data_stream.',1,4);
constellation_matrix=repmat(constellation,length(data_stream),1);

[~,ind] = min(abs(constellation_matrix - demapping_matrix),[],2);
rxbits = de2bi(ind-1, 'left-msb',2);
rxbits = reshape(rxbits.', 1, []);

% If the pseudo-random transformation was active reverse it
if conf.pseudo_activation == 1
    rxbits = bitxor(rxbits, conf.pseudo_sequence);
end

end