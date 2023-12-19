function [txsignal conf] = tx_ofdm(txbits,conf)

% --------------- STAGE A ---------------
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

% --------------- STAGE B ---------------
% Check if divisible by conf.carriers, if not add zeros at the end
if mod(length(mapped_symbols),conf.carriers) ~= 0
    pad = zeros(1,conf.carriers- mod(length(mapped_symbols),conf.carriers));
    mapped_symbols_full = [mapped_symbols pad];
    conf.npad = length(pad);
else
    mapped_symbols_full = mapped_symbols;
end

% --------------- STAGE C ---------------
% Add Training symbol in between OFDM symbols every conf.trainOccurrence [vectorized approach]
if mod(length(mapped_symbols_full)/conf.carriers, conf.trainOccurrence) ~=0
    % Dumb block to allow matrix operation for whatever tranOccurrence value
    n_dumb_blocks = conf.trainOccurrence- mod(length(mapped_symbols_full)/conf.carriers, conf.trainOccurrence);
    dumb_blocks = zeros(1, conf.carriers*n_dumb_blocks);
    mapped_symbols_full = [mapped_symbols_full dumb_blocks];
else
    dumb_blocks = [];
end

mapped_symbols_full = reshape(mapped_symbols_full, conf.carriers*conf.trainOccurrence, []).';


insert_train = ones(size(mapped_symbols_full,1),1).*mapped_training;
matrix_traindata = [insert_train mapped_symbols_full];
matrix_traindata  = reshape(matrix_traindata.', [], 1); % Insert training in proper position
conf.n_training = size(insert_train, 1);

% Remove dumb blocks used to insert training
matrix_traindata = matrix_traindata(1:end-length(dumb_blocks));


% --------------- STAGE D ---------------
% Apply Over-sampled IFFT
%matrix_traindata = osifft(matrix_traindata,conf.os_factor);%*80;
mat_line = []; %zeros(conf.ofdm_nsymbols, conf.carriers*conf.os_factor);
for i=0:length(matrix_traindata)/conf.carriers -1
    mat = osifft(matrix_traindata((i*conf.carriers)+1:(i*conf.carriers)+conf.carriers),conf.os_factor)*80;%*80;
    %mat_line(i+1,:) = mat;% 
    mat_line= [mat_line; mat];
end
matrix_traindata = mat_line;
matrix_traindata = reshape(matrix_traindata, conf.carriers*conf.os_factor, []).';

% matrix_traindata = reshape(matrix_traindata, conf.carriers, []).';
% mat_line = osifft_prova(matrix_traindata, conf.os_factor);



% --------------- STAGE E ---------------
% Cyclic Prefix
CP = matrix_traindata(:, end-conf.Ncp+1:end); % Extract it from the matrix
final_matrix = [CP matrix_traindata]; % Add it

% Apply windowing to each OFDM symbol
% for ii = 1:size(final_matrix, 1)
%     % Apply window to the OFDM symbol
%     windowed_symbol = hamming(conf.carriers*conf.os_factor + conf.Ncp).' .* final_matrix(ii, :);
%     final_matrix(ii, :) = windowed_symbol;
% end

% Add suffix and perform windowing
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


final_vector = reshape(final_matrix.', [], 1).'; % Create final row vector

% --------------- STAGE F ---------------

% Up - Mixing
time = (0:1/conf.f_s:(length(final_vector)-1)/conf.f_s);
conf.time = time;
txsignal = real(final_vector.*exp(2*1j*pi*conf.f_c*time));

conf.len_txsignal = length(txsignal);
%conf.tx_symbols = mapped_symbols_full;


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