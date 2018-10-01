function[est_doa_vec,counts,azg,incg,smoothed_hist] = sspiv_est_doa(z,fs,mic,c,nSrc)

%% deal with optional inputs
if nargin < 4 || isempty(c)
    c = soundspeed; %requires voicebox
end
if nargin < 5 || isempty(nSrc)
    nSrc = 10; %cant imagine trying to find more than that!
end


%% stft
% settings
frame_duration = 0.008;
frame_overlap_frac = 0.25;

% processing
nwin = round(frame_duration * fs);
nfft = nwin;
ninc = round(frame_overlap_frac * frame_duration * fs);
win = hamming(nwin,'periodic');
[Z,pm] = stft(z,win,ninc,nfft,fs);

%% sht + mode strength compensation
N_harm = 3; % spherical harmonic order to use

% analysis range
sht_params.freq_min = 500;
sht_params.freq_max = 4000;

[min_f_bin,max_f_bin] = sh_valid_frequency_bins(mic,pm,N_harm,sht_params); % determine valid frequency range for analysis
pm.f = pm.f([min_f_bin:max_f_bin]);
Praw = sht(Z([min_f_bin:max_f_bin],:,:),mic,N_harm);
Pcomp = compModeStrengthDPD(Praw,pm.f,mic,N_harm,c);

%% sspiv
% time frequency range over which to smooth the spatial covariance
time_span = 0.032;
freq_span = 250;


% processing
% - convert smoothing region to number of bins
tf_smooth_params.nAvgFrames = timespan2frames(fs,nwin,ninc,time_span);
tf_smooth_params.nAvgFreq = freqspan2bins(fs,nfft,freq_span);

[R, R_extras] = tf_smoothed_correlation_v3(Pcomp,pm,mic,tf_smooth_params);
[Ucell,Scell] = svd_for_cells(R);
[sspivs,svr_ratio] = sspiv1(Ucell,Scell);


%% doa est
smooth_sigma = 4;
azres = 2;
incres = 2;

stacked_pivs = reshape(permute(sspivs,[1 3 2]),[],3);
[counts,azg,incg] = piv_hist(stacked_pivs,azres,incres);
[smoothed_hist] = smooth_histogram_pd(counts,smooth_sigma,azg,incg);
[est_doa_vec,~,peak_vals] = sph_local_peaks_to_doa_vecs_with_peak_vals(smoothed_hist,azg,incg,nSrc);

