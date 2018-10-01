% demo_generate_data produces a multichannel wav file that might have 
% been recorded from a particular spherical microphone array
% 
% see README file for external dependencies
% 
% Copyright (C) Alastair H. Moore 2017 

%% Header information
%
% -- Change history --
% Date          Author          Change
% 18/07/2017	A.H. Moore      First version
%
%
% -- SVN version details --
% Last changed in commit: $Rev::               $
%    Date of last commit: $Date::              $
%  Author of last commit: $Author::            $
%
% -- License details --
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
% Directions of arrivals of sources specified in azimuth (az) and
% inclination (inc)
source_az_deg  = [40 120 200].';
source_inc_deg = [80 100 80].';

% Monoaural anechoic speech file associated with each source
source_file_path = fullfile('example_audio',{'F1s3.wav','F3s3.wav','M9s3.wav'});

% Sample rate for processing
fs = 8000;    % Hz

% Duration of output to generate
duration = 5; % seconds

% assumed speed of sound
c = soundspeed;

% Microphone to simulate
% - array is determined by the following fields
%              a: radius of sensor elements (same for all elements)
%   sensor_angle: [nSensors 2] matrix of [azimuth inclination]
%           quad: [nSensors 1] quadrature weight associated with each
%                  sensor for discrete spherical Fourier transform
%        sphType: 'rigid' or 'open' specifies the array configuration
mic = define_eigenmike();

% Length of filter in the time domain (in samples) 
% - determines the number of discrete frequencies at which to evaluate
%   microphone response
filt_len = 4095;

% Spherical harmonic order truncation error
% - representing a plane wave as a truncated series leads to errors in the
%   simulated microphone signals
% - here the maximum allowable error is specified
max_truncation_error_db = -60;

% Output wav file path and word length
out_file_path = 'demo_simulated_recording';
nBits = 24;

%% Processing
% convert to radians for internal processing
source_az_rad = source_az_deg .* pi/180;
source_inc_rad = source_inc_deg .* pi/180;
nSrc = size(source_az_rad,1);

% determine the spherical harmonic order required to obtain good
% approximation to microphone signals
sph_trunc_order = minTruncationOrder(mic,fs,c,max_truncation_error_db);
% or could set it directly to desired value to ensure no spatial aliasing
% sph_trunc_order = 3;

% mode strength of the array configuration
% - determined for each frequency and spherical harmonic order
% - only calculate positive frequencies of fft
nfft = filt_len;
nFreq = fix(nfft/2)+1;
freq_scale = (0:nFreq-1).' * fs/nfft;
%   b: [nFreq,sph_trunc_order+1] - mode stregnth
% i_b: [1 nSH] - indices required to expand (replicate) entries in b to
%      match the number of spherical harmonics
[b,i_b] = modeStrength(mic.sphType,mic.a,mic.a,freq_scale,sph_trunc_order,c);
%b(1,:) = b(2,:); %undefined at dc so use value from first non-zero frequency
b(1,:) = zeros(1,sph_trunc_order+1); %can't have dc gain


% spherical harmonic functions for the source directions
% [nSrc, nSH]
Y_doa = sphBasis(source_az_rad,source_inc_rad,sph_trunc_order);

% spherical harmonic functions for the microphone directions
% [nSensor, nSH]
nSensor = size(mic.sensor_angle,1);
Y_sensor = sphBasis(mic.sensor_angle(:,1),mic.sensor_angle(:,2),sph_trunc_order);

% spherical fourier transform of response to each plane wave
% [nFreq, nSH, nSrc]
Hnm = bsxfun(@times,b(:,i_b),permute(conj(Y_doa),[3 2 1]));

% inverse spherical fourier transform to get freqeuncy response at each
% microphone due to each source
%                         Hnm: [nFreq, nSH, nSrc]
% permute(Y_sensor,[3 2 4 1]): [    1, nSH,    1,   nSensor]
%                           H: [nFreq,   1, nSrc,   nSensor]
H = sum(bsxfun(@times,Hnm,permute(Y_sensor,[3 2 4 1])),2);

% inverse fourier transform to get corresponding impulse responses
% h: [len_filt,  1, nSrc, nSensor]
h = ifft(H,nfft,1,'symmetric');

% rearrange
% h: [len_filt,  nSensor, nSrc]
h = circshift(squeeze(permute(h,[1 4 3 2])),[fix(filt_len/2)-1 0 0]);


% finally load in the speech and convolve with microphone responses
nSamplesRequired = ceil(duration*fs);
mic_signals = zeros(nSamplesRequired,nSensor);

for iSrc = 1:nSrc
    ainfo = audioinfo(source_file_path{iSrc});
    in_fs = ainfo.SampleRate;
    in_sig = audioread(source_file_path{iSrc},[1 ceil(duration*in_fs)]);
    if in_fs~=fs
        in_sig = resample(in_sig,fs,in_fs);
    end
    in_sig = activlev(in_sig,fs,'n');
    tmp_mic_signals = fftfilt(squeeze(h(:,:,iSrc)),in_sig);
    nSamplesAvailable = size(tmp_mic_signals,1);
    if nSamplesAvailable>nSamplesRequired
        mic_signals = mic_signals + tmp_mic_signals(1:nSamplesRequired,:);
    else
        mic_signals(1:nSamplesAvailable,:) = mic_signals(1:nSamplesAvailable,:) + tmp_mic_signals;
    end
end

% add some sensor noise
snr_db = 10;
sig_pow = mean(mic_signals(:).^2);
mic_signals = mic_signals + 10.^(-snr_db/20) * randn(size(mic_signals));

%normalise and write wav file
audiowrite([out_file_path '.wav'],normalise(mic_signals,nBits),fs,'BitsPerSample',nBits)

%save ground truth source directions and microphone array configuration 
save([out_file_path '.mat'],'source_az_deg','source_inc_deg','mic','c');
    


