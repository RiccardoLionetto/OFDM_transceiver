% --------------------------- %
% OFDM Final Project
%
% Informations:
% - Operating mode: 'matlab' (generic MATLAB audio routines)
% --------------------------- %


close all,clear all

%% configuration Values
% Variables needed for:

% ... loading an image to be transmitted
image = imread("image_smallSize.jpeg"); % 'image_fullSize' to use the big version - 'image_smallSize' for small version - 'black_image' for black image
image = rgb2gray(image); % to remove coloring: [W x H x 3] -> [W x H x 1]
[n_rows, n_cols, n_ch] = size(image);
pixel_stream = reshape(image, [n_rows *n_cols *n_ch], 1);
bitstream = reshape(dec2bin(pixel_stream, 8).', 1, []);
bitstream = bitstream - '0'; % Convert binary digits to a single sequence



data_length = length(bitstream) / 2; % Number of QPSK data symbols

% ... OFDM transeiving
conf.carriers = 256;
conf.bitsXsymb = conf.carriers*2; % Because we are using QPSK
conf.spacing = 5; % Spacing between carriers in [Hz]
conf.ofdm_nsymbols = ceil(data_length/conf.carriers); % Number of OFDM tx symbols
conf.f_s     = 48000;   % sampling rate  
conf.os_factor = conf.f_s/(conf.spacing*conf.carriers);
conf.Ncp = 900;%conf.os_factor*conf.carriers/8;
conf.train_length = conf.carriers;
%conf.cpLength = ;
conf.f_c     = 8000;
conf.trainOccurrence = 64; % n° of symbols after which a new training sequence is added
h_hat = 0; % channel estimate (done using training)

% ... to send preamble: single carrier
conf.f_sym   = 100;     % symbol rate
conf.npreamble  = 200;
conf.roll_off = 0.22;
zero_crossing= 5;
conf.filterlength = ceil(conf.os_factor) * zero_crossing;
conf.os_factor_preamble  = conf.f_s/conf.f_sym;

% ... addressing audio on matlab
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;
conf.audiosystem = 'matlab'; % Values: 'matlab','native','bypass'

% ADVANCED CORRECTIONS
% Pseudo random transformation of bitstream
conf.pseudo_activation = 1; % 0 to deactivate
% Suffix and Windowing
conf.suffix = 0; % 1 to activate it
conf.Ncsuffix = 0.05*conf.carriers*conf.os_factor;
% Viterbi-Viterbi
conf.ViterbiViterbi = 1; % 0 to deactivate it
conf.carrierOffset = 0; % To simulate a carrier offset in reception


% Generation of data and preamble and transmission

% DATA ~ QPSK OFDM Symbol
[tx_TrainData conf] = tx_ofdm(bitstream, conf);

% PREAMBLE ~ BPSK Single carrier
preamble = tx_preamble(conf);
%preamble = []; % DEBUG

% Final sequence of transmission (after normalization)
% average_energy_preamble = mean(abs(preamble).^2);
% average_energy_txTrainData = mean(abs(tx_TrainData).^2);
% preamble = preamble/sqrt(average_energy_preamble);
% tx_TrainData = tx_TrainData/sqrt(average_energy_txTrainData);

tx_final = [preamble' tx_TrainData];

tic


for k=1:1
    % ------------------------------------------- %
    % PREAMBLE ~ BPSK Single carrier
%     if k==1
%       preamble = tx_preamble(conf);
%     else
%       preamble = [];   
%     end
%     % ------------------------------------------- %
%     % TRAINING ~ BPSK OFDM Symbol
%     if rem(k,conf.trainOccurrence)==0
%       conf.mapping = 2;
%       txbits_training = randi([0 1],conf.bitsXsymb,1); % Generate
%       [training_txsignal, conf] = tx_ofdm(txbits_training,conf); % Transmit
%     else
%       training_txsignal = [];
%     end
%     % ------------------------------------------- %
%     % DATA ~ QPSK OFDM Symbol
%     conf.mapping = 4;
%     txbits = randi([0 1],conf.bitsXsymb,1); % Generate
%     [txsignal, conf] = tx_ofdm(txbits,conf); % Transmit
% 
% 
%     % Merge (when existing): Preamble + Training + Payload
%     txframe = [preamble;training_txsignal';txsignal'];

    % ------------------------------------------------------------------ %
    % Begin
    % Audio Transmission
    % % % % % % % % % % % %
    
    %rawtxsignal = tx_final; % DEBUG
    % normalize values
%     peakvalue       = max(abs(tx_final));
%     normtxsignal    = tx_final / (peakvalue);
    normtxsignal = tx_final; % DEBUG

    % create vector for transmission
    if k==1 % add initial and final zero padding only when preamble is transmitted
        rawtxsignal = [ zeros(1, conf.f_s) normtxsignal zeros(1 ,conf.f_s)];
    else
        rawtxsignal = normtxsignal;
    end
    % FINO A QUI DECOMMENTARE

    %rawtxsignal = [  rawtxsignal  zeros(size(rawtxsignal)) ]; % add second channel: no signal
    txdur       = length(rawtxsignal)/conf.f_s; % calculate length of transmitted signal
    
%     wavwrite(rawtxsignal,conf.f_s,16,'out.wav')   
    audiowrite('out.wav',rawtxsignal,conf.f_s)  
    
    % Platform native audio mode 
    if strcmp(conf.audiosystem,'native')
        
        % Windows WAV mode 
        if ispc()
            disp('Windows WAV');
            wavplay(rawtxsignal,conf.f_s,'async');
            disp('Recording in Progress');
            rawrxsignal = wavrecord((txdur+1)*conf.f_s,conf.f_s);
            disp('Recording complete')
            rxsignal = rawrxsignal(1:end,1);

        % ALSA WAV mode 
        elseif isunix()
            disp('Linux ALSA');
            cmd = sprintf('arecord -c 2 -r %d -f s16_le  -d %d in.wav &',conf.f_s,ceil(txdur)+1);
            system(cmd); 
            disp('Recording in Progress');
            system('aplay  out.wav')
            pause(2);
            disp('Recording complete')
            rawrxsignal = audioread('in.wav');
            rxsignal    = rawrxsignal(1:end,1);
          end
        
    % MATLAB audio mode
    elseif strcmp(conf.audiosystem,'matlab')
        disp('MATLAB generic');
        playobj = audioplayer(rawtxsignal,conf.f_s,conf.bitsps);
        recobj  = audiorecorder(conf.f_s,conf.bitsps,1);
        record(recobj);
        disp('Recording in Progress');
        playblocking(playobj)
        pause(0.5);
        stop(recobj);
        disp('Recording complete')
        rawrxsignal  = getaudiodata(recobj,'int16');
        rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16'));
        rxsignal = rxsignal.';
        
    elseif strcmp(conf.audiosystem,'bypass')
        rxsignal    = rawtxsignal;
    end
  
    toc
    
    % DA TOGLIERE IN FUTURO DEBUG
    %rxsignal = txframe.';
    %rxsignal = normtxsignal; % skipped zero padding
    % Plot received signal for debugging
    figure;
    plot(rxsignal);
    title('Received Signal')
    
    % End Audio Transmission   
    
    conf.txbits = bitstream;
    [rxbits conf] = rx_ofdm(rxsignal,conf);
    bit_error_rate = sum(rxbits ~= bitstream)/length(bitstream)
    bit_error = sum(rxbits ~= bitstream)

    % ------------------------------------------------------------------ %
% DISPLAY RECEIVED IMAGE
% Convert the bitstream back to a binary matrix
pixel_stream_reconstructed = reshape(rxbits, 8, []).'; % Reshape to 8 bits per pixel
pixel_stream_reconstructed = bin2dec(num2str(pixel_stream_reconstructed)); % Convert binary to decimal

% Reshape the pixel stream back to the original image size
image_reconstructed = reshape(pixel_stream_reconstructed, [n_rows, n_cols, n_ch]);

% Display the reconstructed image
figure;
imshow(uint8(image_reconstructed));
title('Reconstructed Image');

% Compare the original and reconstructed images
figure;
subplot(1, 2, 1);
imshow(image);
title('Original Image');
subplot(1, 2, 2);
imshow(uint8(image_reconstructed));
title('Reconstructed Image');



    % Down mixing
    % rx PREAMBLE
%     if k==1
%       % Call frame synchronizer to get the beginning of data (CP)
%       beginning_of_data = frame_sync(rxsignal, conf);
%       %beginning_of_data = frame_sync(rxsignal.', conf);
%     else
%       beginning_of_data = 1;   
%     end
%     rxsignal = rxsignal(beginning_of_data:end); % cut previous useless part
%     % ------------------------------------------- %
%     % rx TRAINING
%     if rem(k,conf.trainOccurrence)==0
%       % Use the training block for a new estimate of the channel 
%       % Else: h_hat will be the same of the previous frame
%       h_hat = ch_estimation(rxsignal(1:length(training_txsignal)), conf);
%     end
%     % ------------------------------------------- %
%     % rx DATA
%     rxsignal = rxsignal(1+length(training_txsignal):end); % remove training and proceed to demap data
%     [rxbits conf] = rx_ofdm(rxsignal,conf, h_hat);
%     
%     biterrors = sum(rxbits ~= txbits)
% 
%     image = image_decoder(rxbits, image_size);
%     % Compute Bit error rate
%     res.rxnbits(k)      = length(rxbits);  
%     res.biterrors(k)    = sum(rxbits ~= txbits);
    
end