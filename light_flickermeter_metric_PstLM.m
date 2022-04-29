function [PstLM,P_inst_max] = light_flickermeter_metric_PstLM(I, FS)
% 
% Description:
% This function implements a light flickermeter in accordance with IEC TR 61547-1 [1]. 
% The light flickermeter specified in IEC TR 61547-1 [1] can be used to evaluate the flicker severity of light waveforms in terms of the short term flicker metric PstLM.
% This metric represents the flicker perception of light waveforms in a much more objective way than the metrics flicker index or flicker percent (also called modulation depth).
% The function uses some parts of the Matlab function 'flicker_sim' [2] developed by Patrik Jourdan of Solcept [3].
% The function 'flicker_sim' models the voltage fluctuation flickermeter specified in IEC 61000-4-15 [4] to calculate the short-term flicker metric Pst usinga mains voltage waveform as input.
% Further details on the light flickermeter can be found in [5]. General information on Temporal Light Artefacts can be found in [6].
%
% Inputs:
%   I:      vector containg data of the light waveform
%   FS:     sampling frequency of u in Hz (should be >= 2000)
%
% Outputs:
%   PstLM:   short-term flicker metric
%   P_inst_max: max of instantaneous Pst
%
% References:
% [1] IEC TR 61547-1, Equipment for general lighting purposes - EMC immunity requirements - Part 1: An objective voltage fluctuation immunity test method, 2015-04.
% [2] Patrik Jourdan, Flickermeter Simulator, Matlab script 'flicker_sim' posted on the Matlab Central, version 1.6 (Updated 26 Aug 2014).
% [3] Website Solcept: Tools download area: www.solcept.ch/en/tools/flickersim
% [4] IEC 61000-4-15, Electromagnetic compatibility (EMC), Testing and measurement techniques, Flickermeter, Edition 2, 2010-08.
% [5] J.J. Gutierrez, P.A. Beeckman, I. Azcaratea, PROTOCOL TO TEST THE SENSITIVITY OF LIGHTING EQUIPMENT TO VOLTAGE FLUCTUATIONS, CIRED Conference 2015, Lyon, 15-18 June 2015.
% [6] CIE TN 006:2016, CIE Technical note on Visual aspects of time-modulated lighting systems - definitions and measurement models
%
% Notes:
% - This function requires MATLAB with the Signal Processing Toolbox
% - Duration of the light waveform must be > 20 sec
%
%
% COPYRIGHT © PHILIPS LIGHTING HOLDING B.V. 2017
% P. Beeckman 
% Philips Lighting - Eindhoven - Netherlands
% Version: 1.0
% Version: 2017-06-23
%
%% Check input parameters
% check number of input arguments = 2
if (nargin ~= 2)
  error('Invalid number of arguments');
end
% check whether I is a vector
if (~isvector(I))
  error('First input argument must be a vector');
end
% convert vector I to row vector if needed
I = reshape(I, 1, length(I));
% check sample rate
if (FS < 2000)
  warning('Sampling frequency should be >= 2000 Hz');
end
% check whether I has values not below 0
if I<0
   warning('Input light waveform contains negative values! Values must be => 0!')
end
%
%% CALCULATIONS I.A.W. BLOCKS SPECIFIED IN IEC TR 61547-1 [1]
%% Block a: Input illuminance adaptor; see A.2.2 of [1]
u_0 = I / mean(I); % normalises the illuminance by dividing with the average I-level
% 
%% Block b: Filters; see A.2.3 of [1]
% Block 3a: Bandpass filter; equation (A.4) of [1]
% NOTE: this filter is not adapted for different mains frequencies; hence it is applied independent from the mains frequency
% 
% definition of the bandpass filter orders and cut-off frequencies
HIGHPASS_ORDER  = 1;
HIGHPASS_CUTOFF = 0.05;     % high-pass cut-off frequency in Hz
LOWPASS_CUTOFF  = 35;       % low-pass cut-off frequency in Hz
LOWPASS_ORDER   = 6;
%
% subtract DC component to limit filter transients at start of simulation
% (see 'flicker_sim' function [2])
u_0_ac = u_0 - mean(u_0);
%
% application of the band-pass filter
[b_hp, a_hp] = butter(HIGHPASS_ORDER, HIGHPASS_CUTOFF / (FS / 2), 'high');
u_hp = filter(b_hp, a_hp, u_0_ac);

% smooth start of signal to avoid filter transient at start of simulation
smooth_limit = min(round(FS / 10), length(u_hp));
u_hp(1 : smooth_limit) = u_hp(1 : smooth_limit) .* linspace(0, 1, smooth_limit);

[b_bw, a_bw] = butter(LOWPASS_ORDER, LOWPASS_CUTOFF / (FS / 2), 'low');
u_bw = filter(b_bw, a_bw, u_hp);

%% Block 3b: Weighting filter; see equation (A.5) of [1]
% This is the weighting filter of the eye function only; the lamp part has been removed
% Below the parameters of the filter transfer function are given
B1 = 0.041661;
B2 = 44.758;
B3 = 2715.6;
B4 = 29839;
B5 = 0;
A1 = 1;
A2 = 196.32;
A3 = 11781;
A4 = 534820;
A5 = 3505380;
%
num1 = [B1, B2, B3, B4, B5];
den1 = [A1, A2, A3, A4, A5];
[b_w, a_w] = bilinear(num1, den1, FS);
u_w = filter(b_w, a_w, u_bw);
%
%% Block c: Squaring and smoothing; see A.2.4 of [1]; this block is the same as Block 4 of the voltage fluctuation flickermeter [4]
% factor to calibratie P_inst_max to the right value of Table A.1 of IEC TR 61547-1 [1][5]; see equation (A.6) 
% scaling factor derived experimentally by Bilbao University using average of 10 different 60 W inc lamps [5].     
SCALING_FACTOR   =  1.101603155420234e+06;      % scaling of output to perceptibility scale  (according [2])
%
u_q = u_w .^ 2;                                 % squaring of the input signal of block 
%
% definition of the first order low-pass filter (sliding mean filter)
LOWPASS_2_ORDER  = 1;
LOWPASS_2_CUTOFF = 1 / (2 * pi * 300e-3);       % time constant 300 msec
[b_lp, a_lp] = butter(LOWPASS_2_ORDER, LOWPASS_2_CUTOFF / (FS / 2), 'low');
%
s = SCALING_FACTOR * filter(b_lp, a_lp, u_q) ;  % application of the scaling and the low-paas filter
P_inst = s ;                                    % output of block 4 is instantaneous flicker Pinst
%
tau_transient  = 20;                            % transient time in sec that is applied to skip the first part of Pinst to calculate P_inst_max
n_transient    = tau_transient * FS;            % index after the transient from which the max of Pinst is to be calculated
P_inst_max     = max(P_inst(n_transient:end));  % max value of Pinst
%
%% Block d: Statistical evaluation; see A.2.5 of [1]; this block is the same as Block 5 of the voltage fluctuation flickermeter [4]
s=s(n_transient:end);                               % input signal of block d is the instantaneous flicker signal Pinst (transient excluded)
%
NUMOF_CLASSES = 10000;                              % number of bins used for the histogram
[bin_cnt, cpf.magnitude] = hist(s, NUMOF_CLASSES);  % sorts data into the number of bins specified by the scalar NUMOF_CLASSES.
cpf.cum_probability = 100 * (1 - cumsum(bin_cnt) / sum(bin_cnt));
% Calculation of the various percentiles (smoothed)
p_50s = mean([get_percentile(cpf, 30), get_percentile(cpf, 50), get_percentile(cpf, 80)]);
p_10s = mean([get_percentile(cpf, 6), get_percentile(cpf, 8), ...
  get_percentile(cpf, 10), get_percentile(cpf, 13), get_percentile(cpf, 17)]);
p_3s = mean([get_percentile(cpf, 2.2), get_percentile(cpf, 3), get_percentile(cpf, 4)]);
p_1s = mean([get_percentile(cpf, 0.7), get_percentile(cpf, 1), get_percentile(cpf, 1.5)]);
p_0_1 = get_percentile(cpf, 0.1);
%
% Calculation of the final output, i.e. short-term flicker metric PstLM; see 5.7.2 of [4]
PstLM = sqrt(0.0314 * p_0_1 + 0.0525 * p_1s + 0.0657 * p_3s + ...
  0.28 * p_10s + 0.08 * p_50s) ;
%
end  % end of function light_flickermeter_metric_PstLM
%
%% Subfunction: get_percentile
function val = get_percentile(cpf, limit)
  [dummy, idx] = min(abs(cpf.cum_probability - limit));
  val = cpf.magnitude(idx);
end % end of subfunction: get_percentile

