function [after] = ofdmlowpass(before,conf,f)
% LOWPASS lowpass filter
% Low pass filter for extracting the baseband signal 
%
%   before  : Unfiltered signal
%   conf    : Global configuration variable
%   f       : Corner Frequency
%
%   after   : Filtered signal
%
% Note: This filter is very simple but should be decent for most 
% application. For very high symbol rates and/or low carrier frequencies
% it might need tweaking.
%
after = lowpass(before,f,conf.f_s,StopbandAttenuation=30);

