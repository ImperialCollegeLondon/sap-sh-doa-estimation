function[est_doa_vec,srp,az_grid,inc_grid] = est_doas_pwdsrp(z,fs,mic,in_params)
% estimates the directions of arrival of up to nSrcMax sources using a
% steered response power method
% the srp grid is equally spaced (with possibly different resolutions) in 
% azimuth and inclination


%% parameters
% define defaults
params.c = soundspeed;              % speed of sound (m/s)
params.nSrcMax = 10;                % cant imagine trying to find more than that!
params.frame_duration = 0.008;      % stft frame in seconds
params.frame_overlap_frac = 0.75;   % stft frame overlap as fraction
params.az_space_deg = 4;            % angular resolution in azimuth (degrees)
params.inc_space_deg = 4;           % angular resolution in inclination (degrees)
params.sph_harmonic_order_max = []; % the maximum order of SH to compute in the spherical fourier transform
params.freq_min = 500;
params.freq_max = 4000;

% use in_params to override the defaults
if nargin > 3 && ~isempty(in_params)
    params = override_valid_fields(params,in_params);
end

% deal with parameters whose defaults depend on the required inputs
if isempty(params.sph_harmonic_order_max)
    % use some basic logic to infer the maximum order
    % (there are better ways!)
    nSensors = size(mic.sensor_angle,1);
    params.sph_harmonic_order_max = fix(sqrt(nSensors)-1) - 1;
    warning('sph_harmonic_order_max was not specified so guessing...')
    fprintf('Using sph_harmonic_order_max: %d\n',params.sph_harmonic_order_max);
end


%% stft
nwin = round(params.frame_duration * fs);
nfft = nwin;
ninc = round((1-params.frame_overlap_frac) * params.frame_duration * fs);
win = hamming(nwin,'periodic');
[Z,pm] = stft(z,win,ninc,nfft,fs);


%% spherical harmonic transform
% - only do transform for frequencies we want to use

% determine valid frequency range for analysis
% - copy required parameters into a new stucture to pass in
sub_params.c  = params.c;
sub_params.freq_min  = params.freq_min;
sub_params.freq_max  = params.freq_max;
[min_f_bin,max_f_bin] = sh_valid_frequency_bins(mic, pm,...
    params.sph_harmonic_order_max, sub_params);

% remove entries from list of frequencies
pm.f = pm.f(min_f_bin:max_f_bin);

% do the sht, again only using the desired frequencies
Praw = sht(Z(min_f_bin:max_f_bin,:,:),mic,params.sph_harmonic_order_max);


%% mode strength compensation
% - includes frequency depenendent part of beamforming
% - since frequency range of interest is restricted we do not do anything
%   special to control white noise gain

% get the mode strength for mic configuration at each frequency and SH
%   b: [nFreq,sph_trunc_order+1] - mode stregnth
% i_b: [1 nSH] - indices required to expand (replicate) entries in b to
%      match the number of spherical harmonics
[b,i_b] = modeStrength(mic.sphType,mic.a,mic.a,pm.f,...
    params.sph_harmonic_order_max, params.c);
if any(isnan(b(:)))
    error('Need to deal with NaNs in modeStrength')
end

% apply the compensation
%     Praw: [nFreq,nSH,nFrames]
% b(:,i_b): [nFreq,nSH,1]
%    Pcomp: [nFreq,nSH,nFrames]
Pcomp = bsxfun(@rdivide,Praw,b(:,i_b));


%% srp
% - define the grid
% - get steering vectors
% - do beamforming for each direction
% - get power in each direction

% check that spacing leads complete coverage of sphere
if rem(360/params.az_space_deg,1)~=0 || rem(180/params.inc_space_deg,1)~=0
    error('az_space_deg and inc_space_deg should be factors of 360 and 180 respectively')
end
az_vec = pi/180 .* (0:params.az_space_deg:360-params.az_space_deg).';
inc_vec = pi/180 .* (params.inc_space_deg/2 : params.inc_space_deg : 180-params.inc_space_deg/2).';
[inc_grid,az_grid] = ndgrid(inc_vec,az_vec);

% evaluate spherical harmonics in the look directions
% [nDir,nSH]
nDir = numel(inc_grid);
Y_steer = sphBasis(az_grid(:),inc_grid(:),params.sph_harmonic_order_max);

% dimensions of SH domain signal
[nFreq,nSH,nFrames] = size(Pcomp);

% do beamforming and calucalate power
% - freqeuncy dependency already accounted for so we can collapse time and
% frequency dimensions to give a single matrix multiplication
% - limited memory for dense grids makes a single multiplication impossible
%   so a loop is required
% - first get signal into convenient arrangement
% Pcomp: [nSH, nFreq*nFrames]
Pcomp = reshape(permute(Pcomp,[2,1,3]),nSH,nFreq*nFrames);

% - multiplication itself has dimensions
%   [1,nFreq*nFrames] = [1,nSH] * [nSH, nFreq*nFrames]
% - then find mean power over time and frequency
srp_vec = zeros(nDir,1);
reverseStr = '';
tic
for idir = 1:nDir
    if 1 && rem(idir,100)==0
	msg = sprintf('Processing beam %d/%d', idir, nDir);
        fprintf([reverseStr, msg]);
        reverseStr = repmat(sprintf('\b'), 1, length(msg));
    end
    srp_vec(idir) = mean(abs(Y_steer(idir,:) * Pcomp).^2,2);
end
elapsed = toc;
fprintf('...Done in %d seconds!\n', round(elapsed));

% reshape to match the grid
srp = reshape(srp_vec,size(az_grid));

%% doa est
% nSrcMax given as input paramter, or default at top of script, determines
% the maximum number of peaks in the srp that should be returned
% peak_vals could be used to select which of the identified peaks should be
% retained/discarded
[est_doa_vec,~,peak_vals] = sph_local_peaks_to_doa_vecs_with_peak_vals(...
    srp, az_grid, inc_grid, params.nSrcMax);

