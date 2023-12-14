function [rxbits conf] = rx_ofdm(rxsignal,conf)

% Down mixing
time_rx = (0:1/conf.f_s:(length(rxsignal)-1)/conf.f_s);
rx_downmixing = rxsignal.*exp(-2*1j*pi*conf.f_c*time_rx);

figure
plot(abs(fftshift(fft(rx_downmixing))))
title("down-mixed rx signal")

% Frame synchronization to extract the starting index
% [Lowpass filtering is done inside the frame_sync function]
% [An additional lowpass filter will be applied to the OFDM sequence]

idx_start = frame_sync(rx_downmixing,conf);
%idx_start = 1; % DEBUG
rx_OFDM = rx_downmixing(idx_start:end);

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


%Cyclic prefix removal
% OFDM_td_noCP_matrix = OFDM_td_matrix(:,(conf.Ncp*conf.os_factor)+1:end);
OFDM_td_noCP_matrix = OFDM_td_matrix(:,1+conf.Ncp:end);
OFDM_td_noCP_matrix = reshape(OFDM_td_noCP_matrix.', 1, []);
OFDM_td_noCP_matrix = OFDM_td_noCP_matrix.';

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

% Correct them by using channel estimation: every symbol has to be
% corrected by using the channel estimation of the last sent training
% sequence. Some operations are therefore performed first to allow this computation in a matrix form
expanded_ch_matrix = repelem(H_hat, conf.trainOccurrence, 1);

% ------------------------------—------------------------------—
theta_hat = zeros(size(expanded_ch_matrix,1), conf.carriers);   % Estimate of the carrier phase
data_ch_corrected = zeros(size(data_symbols));


if conf.ViterbiViterbi == 1
% Apply viterbi-viterbi algorithm
for train_index = 1:size(H_hat, 1) % Loop over received training symbols
    theta_hat_train = angle(H_hat(train_index,:));
    %for index=1:size(expanded_ch_matrix, 1) % Loop over data
    for index=(train_index-1)*conf.trainOccurrence+1:(train_index-1)*conf.trainOccurrence+conf.trainOccurrence
        if (train_index-1)*conf.trainOccurrence+conf.trainOccurrence > size(expanded_ch_matrix, 1)
            continue
        end
        %theta_hat(index,:) = angle(expanded_ch_matrix(index,:));
        deltaTheta = (1/4*angle(-data_symbols(index,:).^4)).';% + pi/2*(-1:4);
        deltaTheta = repelem(deltaTheta, 1 ,6);
        rotation = repelem((-1:4), conf.carriers,1);
        deltaTheta_rotations = deltaTheta+rotation;

        % Pick values closer to previous estimation (no fast change approx)
        theta_hat_train_rep = repelem(theta_hat_train.', 1,6);
        % Find closest values and corresponding phase
        [~, ind] = min(abs(deltaTheta_rotations-theta_hat_train_rep), [], 2);
        idx = sub2ind(size(deltaTheta_rotations),(1:numel(ind)).',ind);
        theta = deltaTheta_rotations(idx);
    
        % Lowpass filter phase
        if index == (train_index-1)*conf.trainOccurrence+1
            prev_theta = theta_hat_train;
        else
            prev_theta = theta_hat(index-1);
        end
        theta_hat(index,:) = mod(0.01*theta.' + 0.99*prev_theta, 2*pi);
        
        % Phase correction
        data_ch_corrected(index,:) = data_symbols(index,:) .* exp(-1j * theta_hat(index,:));
    end
end
else
data_ch_corrected = data_symbols ./ expanded_ch_matrix;
end
% ------------------------------—------------------------------—


% Code to see how often training symbols should be sent
figure
hold on
for i=1:conf.trainOccurrence
    subplot(4, 5, i);  % Create a 4x5 grid of subplotsdata_ch_corrected
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
plot((0:length(h_hat_power)-1)*(1/conf.f_s), h_hat_power.', 'Marker', 'none', 'LineWidth', 1.5);
title('Power Delay Profile (PDP)');
xlabel('Delay (samples)');
ylabel('Power');


% Demapping QPSK
bistream = conf.txbits; % DEBUG
constellation = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
demapping_matrix = repmat(data_stream.',1,4);
constellation_matrix=repmat(constellation,length(data_stream),1);

[~,ind] = min(abs(constellation_matrix - demapping_matrix),[],2);
rxbits = de2bi(ind-1, 'left-msb',2);
rxbits = reshape(rxbits.', 1, []);

% Eliminate pseudo-random transformation
rxbits = bitxor(rxbits, conf.pseudo_sequence);

end