function[P] = sht(X,mic,N_harm)
%SHT returns the (uncompensated) eigenbeams given the STFT of the
%microphone signals
%
% Implements the spherical harmonic transform for spherical microphone
% array signals as defined in e.g.
% Rafaely2005 eq 6
% Jarret???? book eq 3.7 and 3.9

%TODO: check that the number of harmonics requested is possible for microphone
%specified

% multiply each microphone signal by it's quadrature weight
X = bsxfun(@times, X, mic.quad.');                                         % [nfreq, nchans, nframes]
% evaluate spherical harmonics at sensor angles
Y = conj(sphBasis(mic.sensor_angle(:,1), mic.sensor_angle(:,2), N_harm));  % [nchans, nbasis]
nbasis = size(Y,2);

% perform transform
% equivalent to looping over ii and doing matrix multiply P(:,:,ii) = X(:,:,ii) * Y
[nfreq, nchans, nframes] = size(X);
P = reshape(permute(X,[1 3 2]),nfreq*nframes,nchans) * Y;                  % [nfreq*nframes, nbasis]
P = permute(reshape(P,nfreq,nframes,nbasis),[1 3 2]);                      % [nfreq, nbasis, nframes]
