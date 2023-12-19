function [rxbits conf] = rx_ofdm(rxsignal,conf)


% --------------- STAGE A: down-mixing (+ carrier offset) ---------------
% Down mixing
time_rx = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
rx_downmixing = rxsignal.*exp(-2*1j*pi*(conf.f_c+conf.carrierOffset)*time_rx);
% nella realtà quando si lavora a frequenze più alte (es. 2.4 GHz, il carrier in tx e quello ricevuto in
% rx non avranno la stessa frequenza  (carrier frequency offset). In questo
% caso il viterbi può aiutare -> plot rotazione di fase senza viterbi e con
% viterbi)

figure
plot(abs(fftshift(fft(rx_downmixing))))
title("down-mixed rx signal")


% --------------- STAGE B: frame synchronization ---------------
% Frame synchronization to extract the starting index
% [Lowpass filtering is done inside the frame_sync function]
% [An additional lowpass filter will be applied to the OFDM sequence]
idx_start = frame_sync(rx_downmixing,conf);
%idx_start = 1; % DEBUG
rx_OFDM = rx_downmixing(idx_start:end);

% --------------- STAGE C: lowpass filtering ---------------
% OFDM signal processing
%Lowpass filter
filtered_signal = ofdmlowpass(rx_OFDM,conf,ceil((conf.carriers+1)/2)*conf.spacing);

figure
plot(abs(fftshift(fft(filtered_signal))))
title("down-mixed and filtered rx signal")

%Serial to parallel conversion
ofdm_frame= filtered_signal(1:conf.len_txsignal);
%OFDM_td_matrix = reshape(ofdm_frame, [], conf.carriers*conf.os_factor+conf.Ncp);
OFDM_td_matrix = reshape(ofdm_frame, conf.carriers*conf.os_factor+conf.Ncp, []).';


% --------------- STAGE D: Cyclic Prefix removal ---------------
%Cyclic prefix removal
% OFDM_td_noCP_matrix = OFDM_td_matrix(:,(conf.Ncp*conf.os_factor)+1:end);
OFDM_td_noCP_matrix = OFDM_td_matrix(:,1+conf.Ncp:end);
OFDM_td_noCP_matrix = reshape(OFDM_td_noCP_matrix.', 1, []);
OFDM_td_noCP_matrix = OFDM_td_noCP_matrix.';

% --------------- STAGE E: Over-sampled FFT ---------------
% Frequency domain using FFT
mat_line = [];
for ii=0:length(OFDM_td_noCP_matrix)/(conf.carriers*conf.os_factor) -1
    mat = osfft(OFDM_td_noCP_matrix((ii*conf.carriers*conf.os_factor)+1:(ii*conf.carriers*conf.os_factor)+conf.carriers*conf.os_factor),conf.os_factor);%*80;
    mat_line = [mat_line; mat];
end
OFDM_fd_vector = mat_line;

figure
plot(OFDM_fd_vector, '.')
title("rx signal after osfft [there are data, training at +-1, zero pads]")

OFDM_fd_vector = reshape(OFDM_fd_vector, conf.carriers, []).';


% --------------- STAGE F: Channel estimation ---------------
% CHANNEL ESTIMATION
% extract training symbols [conf.n_training x conf.carriers]
training_symbols = OFDM_fd_vector(1:conf.trainOccurrence+1:end, :); 
% sent training sequence to make the comaprison
sent_training_sequence = (-1 +2*lfsr_training(conf.train_length))';
sent_training_matrix = repmat(sent_training_sequence, [size(training_symbols,1), 1]);
% division to obtain ch estimation
H_hat = training_symbols ./ sent_training_matrix;

% Eliminate training symbols from OFDM_fd_vector to extract data from incoming T+D stream
OFDM_fd_vector(1:conf.trainOccurrence+1:end, :) = [];
data_symbols = OFDM_fd_vector;

figure
plot(data_symbols, '.')
title("rx signal after osfft and removal of training")


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
    final_ch_estimation_matrix(ss,:) = theta_hat_train.*abs(H_hat(train_index,:));
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
    
        % Lowpass filter phase
%         if index == (train_index-1)*conf.trainOccurrence+1
%             prev_theta = theta_hat_train;
%         else
%             prev_theta = theta_hat(index-1, :);
%         end

        % Phase correction (and magnitude from prev. ch symbol)
        if index ~= (train_index-1)*conf.trainOccurrence+1
        theta_hat(index,:) = mod(0.01*theta.' + 0.99*prev_theta, 2*pi);
        data_ch_corrected(index,:) = data_symbols(index,:) .* exp(-1j * theta_hat(index,:)) ./abs(H_hat(train_index,:));
        else
            % The estimation is taken directly from the last training symbol
            theta_hat(index,:) = prev_theta;
            data_ch_corrected(index,:) = data_symbols(index,:) .* exp(-1j * theta_hat(index,:)) ./abs(H_hat(train_index,:));
        end
        % forse basta user prev theta che è già definito in base a quella
        % condizione
        final_ch_estimation_matrix(ss,:) = abs(H_hat(train_index,:)).*exp(1j*theta_hat(index,:));
        %expanded_ch_matrix(ss,:) = abs(H_hat(train_index,:)).*exp(1j*theta_hat(index,:));
        ss = ss+1;
%         if mod(index, 7) == 0
%         % Plot dei primi 256 simboli prima e dopop correzione per
%         % visualizzarla
%         figure
%         hold on
%         plot(data_symbols(index,:), '.')
%         plot(data_ch_corrected(index,:), '.')
%         hold off
%         legend("Symbols before ch correction", "Symbols after Viterbi correction")
%         title(index+ "° group of arrived symbols")
%         end

    end
end

% DEBUG - poi cancellare

% prova_angle = angle(final_ch_estimation_matrix(:,1));
% rad2deg(prova_angle);
% plot(ans)
% 
% 
% tx_symbols_toReshape = conf.tx_symbols;
% tx_symbols = reshape(tx_symbols_toReshape.', conf.carriers, []).';
% diff = rad2deg(angle(tx_symbols(1:8,1)))-rad2deg(angle(data_symbols(1:8,1)))
% rad2deg(theta_hat(1:8,1))
% fine DEBUG
figure
hold on
plot(rad2deg(angle(data_symbols(:,1))))
plot(rad2deg(theta_hat(:,1)))
%plot(1:length(theta_hat(:,1)),rad2deg(angle(tx_symbols(:,1))))
plot(rad2deg(angle(data_ch_corrected(:,1))))
plot(rad2deg(angle(expanded_ch_matrix(:,1))))
hold off
legend('data phase', 'theta correction viterbi', 'phase after correction', 'phase correction no viterbi')
else
data_ch_corrected = data_symbols ./ expanded_ch_matrix;

figure
subplot(1, 2, 1);
plot(angle(expanded_ch_matrix(:,1)))
title("Phase of ch prediction, 1° subcarrier")
subplot(1, 2, 2);
plot(angle(data_symbols(:,1)))
title("Phase of arrived data symbols, 1° subcarrier")

figure
hold on
plot(angle(expanded_ch_matrix(:,1)))
plot(angle(data_symbols(:,1)))
hold off
legend('ch prediction', 'data symbol')

end
% ------------------------------—------------------------------—


% Code to see how often training symbols should be sent
figure
hold on
n_rows_plot = ceil(conf.trainOccurrence/4);
for i=1:conf.trainOccurrence
    subplot(n_rows_plot, 4, i);  % Create a 4x5 grid of subplotsdata_ch_corrected
    plot(data_ch_corrected(i,:), '.');
    ylim([-2 +2]);
    xlim([-2 +2]);
    %legend(i)
end
hold off

data_stream = reshape(data_ch_corrected.', 1, []);
data_stream = data_stream(1:end-conf.npad);% remove pad added in transmission to have a even n° of symbols

figure
plot(data_stream, '.')
title("rx signal after:osfft, removal of training & zero pads, ch correction")

%Power delay profile
h_hat_fft = fft(H_hat(1,:));
h_hat_power = abs(h_hat_fft).^2;
figure;
%plot(1:length(h_hat_power), h_hat_power);
%plot((0:length(h_hat_power)-1)*(1/conf.f_s), h_hat_power);
plot((0:1/conf.f_s:(length(h_hat_power)-1)/conf.f_s), h_hat_power.', 'Marker', 'none', 'LineWidth', 1.5);
title('Power Delay Profile (PDP)');
xlabel('Delay [second]');
ylabel('Power');

%(0:1/conf.f_s:(length(h_hat_power)-1)/conf.f_s);

% --------------- STAGE G: Demapping ---------------

% Demapping QPSK
bistream = conf.txbits; % DEBUG
constellation = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
demapping_matrix = repmat(data_stream.',1,4);
constellation_matrix=repmat(constellation,length(data_stream),1);

[~,ind] = min(abs(constellation_matrix - demapping_matrix),[],2);
rxbits = de2bi(ind-1, 'left-msb',2);
rxbits = reshape(rxbits.', 1, []);

if conf.pseudo_activation == 1
% Eliminate pseudo-random transformation
    rxbits = bitxor(rxbits, conf.pseudo_sequence);
end


end