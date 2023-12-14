function [txsignal conf] = tx_ofdm(txbits,conf)

% TRAINING
mapped_training = -1 +2*lfsr_training(conf.carriers).'; % BPSK [1xconf.carriers]

% DATA
% In transmission there are peaks in the signal, because some sequences have all equal bits,
% therefore obtaining a peak in the ifft -> apply pseudorandom transformation
if conf.pseudo_activation == 1
    [txbits_ps, conf] = pseudo_rand(txbits, conf);
else
    txbits_ps = txbits;
end
mapped_symbols = QPSK(txbits_ps); % QPSK [dimensions 1x ~]

% Check if divisible by conf.carriers, if not add zeros at the end
if mod(length(mapped_symbols),conf.carriers) ~= 0
    pad = zeros(1,conf.carriers- mod(length(mapped_symbols),conf.carriers));
    mapped_symbols_full = [mapped_symbols pad];
    conf.npad = length(pad);
end

% Add Training symbol in between OFDM symbols every conf.trainOccurrence [vectorized approach]
if mod(length(mapped_symbols_full)/conf.carriers, conf.trainOccurrence) ~=0
    % Dumb block to allow matrix operation for whatever tranOccurrence value
    n_dumb_blocks = conf.trainOccurrence- mod(length(mapped_symbols_full)/conf.carriers, conf.trainOccurrence);
    dumb_blocks = zeros(1, conf.carriers*n_dumb_blocks);
    mapped_symbols_full = [mapped_symbols_full dumb_blocks];
else
    dumb_blocks = 0;
end

mapped_symbols_full = reshape(mapped_symbols_full, conf.carriers*conf.trainOccurrence, []).';


insert_train = ones(size(mapped_symbols_full,1),1).*mapped_training;
matrix_traindata = [insert_train mapped_symbols_full];
matrix_traindata  = reshape(matrix_traindata.', [], 1); % Insert training in proper position
conf.n_training = size(insert_train, 1);

% Remove dumb blocks used to insert training
matrix_traindata = matrix_traindata(1:end-dumb_blocks);

% Apply Over-sampled IFFT
%matrix_traindata = osifft(matrix_traindata,conf.os_factor);%*80;
mat_line = [];
for i=0:length(matrix_traindata)/conf.carriers -1
    mat = osifft(matrix_traindata((i*conf.carriers)+1:(i*conf.carriers)+conf.carriers),conf.os_factor)*80;
    mat_line = [mat_line; mat];
end

matrix_traindata = mat_line;
matrix_traindata = reshape(matrix_traindata, conf.carriers*conf.os_factor, []).';


% Cyclic Prefix
CP = matrix_traindata(:, end-conf.Ncp+1:end); % Extract it form the matrix
final_matrix = [CP matrix_traindata]; % Add it

% Apply windowing to each OFDM symbol
% for ii = 1:size(final_matrix, 1)
%     % Apply window to the OFDM symbol
%     windowed_symbol = hamming(conf.carriers*conf.os_factor + conf.Ncp).' .* final_matrix(ii, :);
%     final_matrix(ii, :) = windowed_symbol;
% end

% Add suffix and perform windowing: scritta da me
if conf.suffix == 1
    final_portion = final_matrix(1, end-conf.Ncsuffix+1:end);
    hamming_1to0 = hamming(2*conf.Ncsuffix);
    hamming_0to1 = hamming(2*conf.Ncp);
    sinc_0to1 = sinc(-pi:1/conf.Ncp:0);
    sinc_1to0 = sinc(0:1/conf.Ncsuffix:pi);
%    final_portion_rrc = [final_portion.*sinc_1to0(1:conf.Ncsuffix) zeros(1, conf.Ncp-conf.Ncsuffix)]; % rrc goes from 1 to 0
    final_portion_rrc = [final_portion.*hamming_1to0(end-conf.Ncsuffix+1:end)' zeros(1, conf.Ncp-conf.Ncsuffix)]; % rrc goes from 1 to 0
    for ii=2:size(final_matrix, 1)
        % Get CP and apply rrc [transition from 0 to 1]
        %CP_rrc = final_matrix(ii, 1:conf.Ncp).*sinc_0to1(end-conf.Ncp+1:end); % rrc goes from 0 to 1
        CP_rrc = final_matrix(ii, 1:conf.Ncp).*hamming_0to1(1:conf.Ncp)';
        windowing = final_portion_rrc+CP_rrc;
        final_matrix(ii, 1:conf.Ncp) = windowing; % Substitute old CP with the new windowed one

        % Get suffix to be used on the next CP and apply rrx [transition from 1 to 0]
        final_portion = final_matrix(ii, end-conf.Ncsuffix+1:end);
        %final_portion_rrc = [final_portion.*sinc_1to0(1:conf.Ncsuffix) zeros(1, conf.Ncp-conf.Ncsuffix)];
        final_portion_rrc = [final_portion.*hamming_1to0(end-conf.Ncsuffix+1:end)' zeros(1, conf.Ncp-conf.Ncsuffix)];
    end
end

% Add suffix and perform windowing: scritta da CHAT
% if conf.suffix == 1
%     final_portion = final_matrix(1, end-conf.Ncsuffix+1:end);
%     hamming_window = hamming(conf.Ncsuffix + conf.Ncp);
%     final_portion_windowed = final_portion .* hamming_window(1:conf.Ncsuffix);
% 
%     for ii = 2:size(final_matrix, 1)
%         % Get CP and apply hamming window [transition from 0 to 1]
%         CP_windowed = final_matrix(ii, 1:conf.Ncp) .* hamming_window(end-conf.Ncp+1:end);
%         windowed_symbol = [CP_windowed final_matrix(ii, conf.Ncp+1:end)];
% 
%         % Substitute old CP with the new windowed one
%         final_matrix(ii, :) = windowed_symbol;
% 
%         % Get suffix to be used on the next CP and apply hamming window [transition from 1 to 0]
%         final_portion = final_matrix(ii, end-conf.Ncsuffix+1:end);
%         final_portion_windowed = final_portion .* hamming_window(1:conf.Ncsuffix);
%     end
% end


final_vector = reshape(final_matrix.', [], 1).'; % Create final row vector

% Up - Mixing
time = (0:1/conf.f_s:(length(final_vector)-1)/conf.f_s);
conf.time = time;
txsignal = real(final_vector.*exp(2*1j*pi*conf.f_c*time));

conf.len_txsignal = length(txsignal);


% Matrix formation: add training in between data every conf.trainOccurrence
% blocks of conf.carriers of data symbols
% TrainData_seq = [];
% num_repeats = floor(length(mapped_symbols_full) / (conf.carriers*conf.trainOccurrence));
% for i=1:num_repeats
%     start_idx = (i-1)*conf.trainOccurrence*conf.carriers+1;
%     data_sym = mapped_symbols_full(start_idx:start_idx+conf.carriers*conf.trainOccurrence-1);
%     TrainData_seq = [TrainData_seq mapped_training data_sym];
% end


% Initialize the result vector
% result_vector = zeros(size(mapped_symbols_full,2) + length(mapped_training) * num_repeats, 1);
% 
% result_vector(1:length(mapped_training)) = mapped_training; % Insert training sequence at the beginning
% % Insert the training sequence every conf.carriers*conf.trainOccurrence symbols
% for i = 1:num_repeats
%     start_idx = length(mapped_training) * i + 1;
%     end_idx = start_idx + length(mapped_training) - 1;
%     result_vector(start_idx:end_idx) = mapped_training;
% end
% % Insert the original symbols
% original_start_idx = length(mapped_training) * (num_repeats + 1) + 1;
% result_vector(original_start_idx:end) = mapped_symbols_full;
% 
% 
% 
% matrix_TrainData = reshape(mapped_symbols_256, conf.trainOccurrence, []);
% matrix_TrainData1(conf.trainOccurrence+1, :) = mapped_training;
% matrix_TrainData2 = reshape(matrix_TrainData1, 1, []);
% matrix_TrainData = [mapped_training.' matrix_TrainData2];
% 
% matrix_TrainData = reshape(matrix_TrainData, [], conf.carriers);
% 
% 
% 
% 
% 
% 
% % lunghezza lowpass filtro ceil((conf.carriers+1)/2)*conf.spacing
% 
% 
% 
% if conf.mapping == 4 % this will be called when generating data
%     % Map bits in QPSK
%     symbols = QPSK(txbits); % Vector of dimensions 256x1 in complex domain
% elseif conf.mapping == 2 % this will be called when generating training
%     % Map bits in BPSK
%     symbols = BPSK(txbits);
% else
%     disp('ERROR: mapping is possible only on BPSK and QPSK'); 
% end
% 
% % Over-sampled IFFT
% OFDM_symbols = osifft(symbols,conf.os_factor)*2;
% 
% % Add Cyclic Prefix
% CP_OFDM_symb = zeros(1, length(OFDM_symbols)+conf.Ncp);
% CP = OFDM_symbols(end-conf.Ncp+1: end);
% CP_OFDM_symb(1:conf.Ncp) = CP;
% CP_OFDM_symb(conf.Ncp+1:end) = OFDM_symbols;
% 
% % Mixing
% time = (0:1/conf.f_s:(length(CP_OFDM_symb)-1)/conf.f_s);
% conf.time = time;
% txsignal = real(CP_OFDM_symb).*cos(2*pi*conf.f_c.*time)-imag(CP_OFDM_symb).*sin(2*pi*conf.f_c.*time);
% 
% % Normalizing signal with respect to average energy
% average_energy = norm(txsignal)^2 /length(txsignal);
% txsignal = txsignal/sqrt(average_energy);
% 
% figure
% plot(time, txsignal)
% title("OFDM datablock sent signal in time")

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
    txbits = reshape(txbits,2,[])';
    GrayMap = 1/sqrt(2) * [(-1-1j) (-1+1j) ( 1-1j) ( 1+1j)];
    symbol = GrayMap(bi2de(txbits, 'left-msb')+1);
end

function[symbol] = BPSK (txbits)
    symbol = -1 + 2*txbits;
end

function [txbits_ps conf] = pseudo_rand(txbits, conf)
conf.pseudo_sequence = randi([0, 1], size(txbits));
txbits_ps = bitxor(txbits, conf.pseudo_sequence);
end

end