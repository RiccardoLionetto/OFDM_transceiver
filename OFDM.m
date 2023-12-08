% --------------------------- %
% OFDM Final Project
%
% Informations:
% - Operating mode: 'matlab' (generic MATLAB audio routines)
% --------------------------- %

% Things to be checked/corrected
printt = "Guardare figura recevied signal dopo wav, ampiezza preamble>>dati -> normalizzare"

% Git https: https://github.com/RiccardoLionetto/ofdm.git
% GIt ssh: git@github.com:RiccardoLionetto/ofdm.git

close all,clear all

%% configuration Values
% Variables needed for:

% ... OFDM transeiving
conf.carriers = 256;
conf.bitsXsymb = conf.carriers*2; % Because we are using QPSK
conf.spacing = 5; % Spacing between carriers in [Hz]
conf.ofdm_nsymbols = 1; % Number of OFDM tx symbols
conf.f_s     = 48000;   % sampling rate  
conf.os_factor = conf.f_s/(conf.spacing*conf.carriers);
conf.Ncp = conf.os_factor*conf.carriers/2;
conf.train_length = conf.carriers;
%conf.cpLength = ;
conf.f_c     = 8000;
conf.trainOccurrence = 7; % n° of symbols after which a new training sequence is added
h_hat = 0; % channel estimate (done using training)

% ... to send preamble: single carrier
conf.f_sym   = 100;     % symbol rate
conf.npreamble  = 100;
conf.roll_off = 0.22;
zero_crossing= 5;
conf.filterlength = ceil(conf.os_factor) * zero_crossing;
conf.os_factor_preamble  = conf.f_s/conf.f_sym;

% ... addressing audio on matlab
conf.bitsps     = 16;   % bits per audio sample
conf.offset     = 0;
conf.audiosystem = 'bypass'; % Values: 'matlab','native','bypass'


for k=1:conf.ofdm_nsymbols+1
    % ------------------------------------------- %
    % PREAMBLE ~ BPSK Single carrier
    if k==1
      preamble = tx_preamble(conf);
    else
      preamble = [];   
    end
    % ------------------------------------------- %
    % TRAINING ~ BPSK OFDM Symbol
    if rem(k,conf.trainOccurrence)==0
      conf.mapping = 2;
      txbits_training = randi([0 1],conf.bitsXsymb,1); % Generate
      [training_txsignal, conf] = tx_ofdm(txbits_training,conf); % Transmit
    else
      training_txsignal = [];
    end
    % ------------------------------------------- %
    % DATA ~ QPSK OFDM Symbol
    conf.mapping = 4;
    txbits = randi([0 1],conf.bitsXsymb,1); % Generate
    [txsignal, conf] = tx_ofdm(txbits,conf); % Transmit


    % Merge (when existing): Preamble + Training + Payload
    txframe = [preamble;training_txsignal';txsignal'];

    % ------------------------------------------------------------------ %
    % Begin
    % Audio Transmission
    % % % % % % % % % % % %
    
    % normalize values
    peakvalue       = max(abs(txframe));
    normtxsignal    = txframe / (peakvalue + 0.3);
    
    % create vector for transmission
    if k==1 % add initial zero padding only when preamble is transmitted
        rawtxsignal = [ zeros(conf.f_s,1) ; normtxsignal];
    else
        rawtxsignal = normtxsignal;
    end

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
        rxsignal     = double(rawrxsignal(1:end))/double(intmax('int16')) ;
        
    elseif strcmp(conf.audiosystem,'bypass')
        rawrxsignal = rawtxsignal(:,1);
        rxsignal    = rawrxsignal;
    end
    
    
    % DA TOGLIERE IN FUTURO
    rxsignal = txframe.';
    %rxsignal = normtxsignal; % skipped zero padding
    % Plot received signal for debugging
    figure;
    plot(rxsignal);
    title('Received Signal')
    
    % End Audio Transmission   
    
    % Now the signal with zero padding is arriving, do we need it?
    % Especially the one at the end...unuseful?
    % ------------------------------------------------------------------ %
   
    % Down mixing
    % rx PREAMBLE
    if k==1
      % Call frame synchronizer to get the beginning of data (CP)
      beginning_of_data = frame_sync(rxsignal.', conf);
      %beginning_of_data = frame_sync(rxsignal.', conf);
    else
      beginning_of_data = 1;   
    end
    rxsignal = rxsignal(beginning_of_data:end); % cut previous useless part
    % ------------------------------------------- %
    % rx TRAINING
    if rem(k,conf.trainOccurrence)==0
      % Use the training block for a new estimate of the channel 
      % Else: h_hat will be the same of the previous frame
      h_hat = ch_estimation(rxsignal(1:length(training_txsignal)), conf);
    end
    % ------------------------------------------- %
    % rx DATA
    rxsignal = rxsignal(1+length(training_txsignal):end); % remove training and proceed to demap data
    [rxbits conf] = rx_ofdm(rxsignal,conf, h_hat);
    
    biterrors = sum(rxbits ~= txbits)

    % Compute Bit error rate
    res.rxnbits(k)      = length(rxbits);  
    res.biterrors(k)    = sum(rxbits ~= txbits);
    
end