% --------------------------- %
% OFDM Final Project
%
% Implementation of an OFDM transceiver
% --------------------------- %


close all,clear all

%% configuration Values
% Variables needed for:

% ... loading an image to be transmitted:
% Three possibilites based on the images on the folder:
% 1. 'image_smallSize.jpeg': standard choice, normal time required to tx [60x68x1 takes around 20 seconds to be transmitted]
% 2. 'image_fullSize.jpeg': higher resolution version of the above image. Extremely long time required to tx
% 3. 'black_image.jpg': all black image to show PAPR problem
image = imread("image_smallSize.jpeg");
image = rgb2gray(image); % to remove coloring: [W x H x 3] -> [W x H x 1]
[n_rows, n_cols, n_ch] = size(image);
pixel_stream = reshape(image, [n_rows *n_cols *n_ch], 1);
bitstream = reshape(dec2bin(pixel_stream, 8).', 1, []);
bitstream = bitstream - '0'; % Convert binary digits to a single sequence

% Uncomment to transmit a random sequence of bits instead of an image
%bitstream = randi([0 1], 1, 500000);

% ... OFDM transceiver
conf.f_s     = 48000;   % sampling rate  
conf.f_c     = 8000; % Carrier frequency
    % OFDM Data symbols
data_length = length(bitstream) / 2; % Number of QPSK data symbols
conf.carriers = 256; % Number of sub-carriers
conf.bitsXsymb = conf.carriers*2; % Because we are using QPSK % DEBUG
conf.spacing = 5; % Spacing between carriers in [Hz]
conf.ofdm_nsymbols = ceil(data_length/conf.carriers); % Number of OFDM symbols transmitted
conf.os_factor = conf.f_s/(conf.spacing*conf.carriers); % Used on the over-sampled IFFT and FFT
conf.Ncp = conf.os_factor*conf.carriers/8; % Length of Cyclic Prefix
    % OFDM Training symbols
conf.train_length = conf.carriers; % Number of sub-carriers of the traning symbols
conf.trainOccurrence = 4; % n° of Data symbols after which a new training sequence is added
h_hat = 0; % channel estimate (done using training)

% ... to send preamble: single carrier
conf.f_sym   = 100;     % symbol rate
conf.npreamble  = 100; % Preamble length (bits)
conf.roll_off = 0.22;
zero_crossing= 5;
conf.filterlength = ceil(conf.os_factor) * zero_crossing;
conf.os_factor_preamble  = conf.f_s/conf.f_sym;

% ... addressing audio on matlab
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;
conf.audiosystem = 'matlab'; % Values: 'matlab','native','bypass'

% ----------- ADVANCED CORRECTIONS -------------
% These can be activated setting thier value to 1.

% - Pseudo random transformation of bitstream: helps reducing the peaks in the 
%   tx signal coming form the IFFT
conf.pseudo_activation = 1; % 1/0 to deactivate
% Suffix and Windowing
conf.suffix = 0; % 1 to activate it
conf.Ncsuffix = 0.05*conf.carriers*conf.os_factor;
% Viterbi-Viterbi
conf.ViterbiViterbi = 0; % 0 to deactivate it
conf.carrierOffset = 0; % To simulate a carrier offset in reception


% --------- Generation of data and preamble and transmission --------------

% DATA ~ QPSK OFDM Symbol
[tx_TrainData conf] = tx_ofdm(bitstream, conf);

% PREAMBLE ~ BPSK Single carrier
preamble = tx_preamble(conf)*3;

tx_final = [preamble' tx_TrainData];

% ------------------------------------------------------------------ %
% Begin Audio Transmission
tic

% Add zero padding at the extremes:
% initial padding: to leave marging from beginning of tx and beginning
% of preamble -> frame synchronization is tested
% final padding: to cover mismatch created by eventual wrong synchronization
rawtxsignal = [ zeros(1, conf.f_s) tx_final zeros(1 ,conf.f_s)];


audiowrite('out.wav',rawtxsignal,conf.f_s)  
% MATLAB audio mode
if strcmp(conf.audiosystem,'matlab')
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
% ----------------- End Audio Transmission ----------------------

% ---------------------- Reception ---------------------------
% Plot received signal for debugging
figure;
plot(rxsignal);
xlabel('Samples')
ylabel('Amplitude')
title('Received Signal - time domain')

% Receiver function
[rxbits conf] = rx_ofdm(rxsignal,conf);


% ----------------------- RESULTS ---------------------------
% Bit error rate and overall number of wrong bits 
bit_error_rate = sum(rxbits ~= bitstream)/length(bitstream)
bit_error = sum(rxbits ~= bitstream)


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