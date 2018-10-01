function[R,extras] = tf_smoothed_covariance(X,pm,nAvgT,nAvgF)
% TF_SMOOTHED_COVARIANCE takes STFT domain signal and computes the outer 
% product over the "channels" averaged across local time and frequency region
% "channels" can be space domain (i.e. microphone signals) or spherical
% harmonic domain (i.e. eigenbeams)

%TODO: Add support for hops in time and frequency

%%
nAvgTot = nAvgT * nAvgF;
[nFreq,nChans,nFrames] = size(X);


%% processing
% two step process
% - find outer product of channels for each TF-bin
% - average over required TF-regions

% outer product
op = reshape(permute(X,[2 4 1 3]), nChans,1,nFreq*nFrames);
op = bsxfun(@times,op, conj(permute(op,[2 1 3])));  %[nChans, nChans, nFreq*nFrames]
op = reshape(op,[nChans, nChans, nFreq,nFrames]);

%preallocate before loop
R = cell(nFreq,nFrames);
extras.valid_bins = zeros(nFreq,nFrames); % gets a one when the corresponding entry of R is filled
extras.t = zeros(nFrames,1); % holds the time at the middle of the corresponding entry of R
extras.f = zeros(nFreq,1);   % holds the frequency at the middle of the corresponding entry of R

%putting result in cell array for historical reasons
for i_frame = 1:(nFrames-nAvgT+1)
    frame_idc = ((i_frame)+[0:nAvgT-1]);
    extras.t(i_frame) = mean(pm.t(frame_idc)); %inefficient but clear
    for i_freq = 1:(nFreq-nAvgF+1)
        freq_idc = ((i_freq)+[0:nAvgF-1]);
        R{i_freq,i_frame} = 1/nAvgTot * sum(sum(op(:,:,freq_idc,frame_idc),3),4);
        extras.valid_bins(i_freq,i_frame) = 1;
        if i_frame == 1
            extras.f(i_freq) = mean(pm.f(freq_idc));
        end
    end
end
