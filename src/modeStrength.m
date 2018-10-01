function[b,i_b] = modeStrength(sphType,a,r,freq,N_harm,c)

%   Mode strengths for plane waves
%
%   [b,i_b] = mode_strength(sphType,a,r,f,N_harm,c)
%
%   Inputs:
%       sphType     Type of sphere ('rigid' or 'open') (default='rigid').
%       a           Radius of sphere (m).
%       r           Evaluation radius (m).
%       freq        Frequencies at which to calculate mode strength (Hz)
%       N_harm      Spherical harmonics order
%       c           Speed of sound (default=343) m/s
%
%   Outputs:
%       b           length(freq) x N_harm+1 matrix of mode strengths
%       i_b         row vector of indices for expanding b
%
%   References:
%       B. Rafaely, "Analysis and Design of Spherical Microphone Arrays,"
%       IEEE Trans. Speech and Audio Processing, vol. 13, no. 1, Jan. 2005.
%
%**************************************************************************
% Author:           M. R. P. Thomas, E. A. P. Habets and D. Jarrett
% Date:             7 November 2012
% Version: $Id: sphModeStr.m 2857 2013-03-29 11:54:42Z dpj05 $
%**************************************************************************

narginchk(5,6);

if(nargin<6)
    c=343;
end
freq = freq(:).'; %force to be row vector internally
nfreqs = length(freq);
b  = zeros(N_harm+1,nfreqs);
farfield_mode_strength = zeros(N_harm+1,nfreqs);

kk = 1 : nfreqs;
lambda = c./freq;
k = 2*pi./lambda;
l = 0 : N_harm;

if strcmp('rigid', sphType)   % rigid sphere
    besselh_derivative = repmat(sqrt(pi./(2*k*a)),N_harm+1,1) .* ( (bsxfun(@besselh,l+0.5-1,k'*a).' - bsxfun(@besselh,l+0.5+1,k'*a).')/2 - bsxfun(@besselh,l+0.5,k'*a).' ./ repmat((2*k*a),N_harm+1,1) );
    if (r == a)
        farfield_mode_strength = 1i./(besselh_derivative .* repmat((k*a).^2,N_harm+1,1));
    else
        besselj_derivative = repmat(sqrt(pi./(2*k*a)),N_harm+1,1) .* ( (bsxfun(@besselj,l+0.5-1,k'*a).' - bsxfun(@besselj,l+0.5+1,k'*a).')/2 - bsxfun(@besselj,l+0.5,k'*a).' ./ repmat((2*k*a),N_harm+1,1) );
        farfield_mode_strength = repmat(sqrt(pi./(2*k*r)),N_harm+1,1) .* bsxfun(@besselj,l+0.5,k'*r).' - besselj_derivative ./ besselh_derivative .* repmat(sqrt(pi./(2*k*r)),N_harm+1,1) .* bsxfun(@besselh,l+0.5,k'*r).';
    end
else                             % open sphere
    farfield_mode_strength = repmat(sqrt(pi./(2*k*r)),N_harm+1,1) .* bsxfun(@besselj,l+0.5,k'*r).';
end
% By including conj() below, we use spherical Hankel functions of the SECOND kind
b = repmat((4*pi*1i.^l).',1,nfreqs) .* conj(farfield_mode_strength);

b = b.'; % makes more sense for compatability with stft output
i_b = [];
for l=0:N_harm
    i_b = [i_b, l+1*ones(1,length(-l:l))];
end