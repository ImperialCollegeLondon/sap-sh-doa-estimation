function[min_f_bin,max_f_bin] = sh_valid_frequency_bins(mic,pm,N_harm_max,in_params)

params.freq_min = [];
params.freq_max = [];
params.c = soundspeed;

params = override_valid_fields(params,in_params);

if isempty(params.freq_min) || isempty(params.freq_max)
    error('params structure must specify freq_min and freq_max')
end

% define the frequency range to use for analysis
% - upper limit due to microphone array
kr = (2 * pi *mic.a / params.c) * pm.f; %wavenumber k = omega/c = 2pi*f/c; r = mic radius
max_f_bin_sh = sum(kr<N_harm_max);                 %nice way of finding last entry in kr to fit below N_harm_max
if max_f_bin_sh == length(pm.f)
    fprintf('All frequency bins are below maximum allowed by kr limit')
end 

% - upper limit imposed by user
max_f_bin = sum(pm.f < params.freq_max);                 %nice way of finding last entry in f to fit below freq_max
if max_f_bin == length(pm.f)
    fprintf('All frequency bins are below maximum allowed by freq_max limit')
end

max_f_bin = min(max_f_bin_sh,max_f_bin);

% - lower limit imposed by user
min_f_bin = sum(pm.f < params.freq_min)+1;    %find first bin to exceed freq_min
if min_f_bin >= max_f_bin
    error('Frequency range is zero!')
end