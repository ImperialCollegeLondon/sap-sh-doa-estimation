% demo_sspiv loads a spherical microphone array recording and uses 
% subspace pseudointensity vectors to estimate the directions of arrival
%
% requires data saved in the form demonstrated in demo_generate_data.m
% 
% see README file for external dependencies
% 
% Copyright (C) Alastair H. Moore 2017 

%% File information
%
% -- Change history --
% Date          Author(s)       Change
% 19/07/2017	A.H. Moore      First version
%
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
% data to use
in_file_path = 'demo_simulated_recording';

% choose the maximum spherical harmonic order to use in analysis
% - this should be chosen according to the microphone array and the maximum
%   frequency of interest
sph_harmonic_order_max = 3;


%% Processing
% load in ground truth data:
%  - source_az_deg
%  - source_inc_deg
% and necessary data required for processing of signals
%  - mic
%  - c
load([in_file_path '.mat'],'source_az_deg','source_inc_deg','mic','c')

% read in the audio data
[z,fs] = audioread([in_file_path '.wav']);

% use plane wave decomposition steered response power to estimate DOAs
% - optional parameter values are passed in as a structure
params.c = c; 
params.sph_harmonic_order_max = sph_harmonic_order_max;

% est_doa_vec: [M 3] Up to M source DOAs represented as cartesian vectors
%         srp: [nInc nAz] power map (summed over frequency range set in function)
%         azg: [nInc nAz] grid of azimuth values
%        incg: [nInc nAz] grid of inclination values
tic;
[est_doa_vec,counts,azg,incg,smoothed_histogram] = est_doas_sspiv(z,fs,mic,params);
elapsed = toc;
fprintf('SSPIV processing done in %2.1f seconds\n',elapsed);

% convert ground truth DOAs into vector form
nSrc = size(source_az_deg,1);
gt_doa_vec = zeros(nSrc,3);
[gt_doa_vec(:,1), gt_doa_vec(:,2), gt_doa_vec(:,3)] = ...
    mysph2cart(pi/180 * source_az_deg, pi/180 * source_inc_deg, ones(nSrc,1));

% plot the estimated results along with ground truth
contour_levels = []; % specify list of levels here to get contour map instead
figure;
map = 10*log10(smoothed_histogram./max(smoothed_histogram(:)));
[~, cax] = plot_2D_map_with_est_gt_doa(azg,incg,map,est_doa_vec,gt_doa_vec,contour_levels);
ylabel(cax,'Normalised smoothed histogram counts [dB]')
set(gca,'clim',[-12 0])
title('Smoothed histogram of SSPIVs with estimated and true DOAs')



