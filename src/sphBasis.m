function Y = sphBasis(az, inc, sphOrdMax)

%   Basis functions for spherical harmonics decomposition evaluated at
%   specified points
%
%   B = sphBasis(azimuth, inclination, N_harm)
%
%   Inputs:
%          az   Mx1 vector of azimuth angles in [0 2pi).
%         inc   Mx1 vector of inclination angles in [0 pi].
%   sphOrdMax   Maximum spherical harmonic order.
%
%   Outputs:
%           B   MxsphOrdMax^2+2*sphOrdMax+1 matrix of spherical harmonics
%               evaluated at (inc, az).
%
%   References:
%       B. Rafaely, "Analysis and Design of Spherical Microphone Arrays,"
%       IEEE Trans. Speech and Audio Processing, vol. 13, no. 1, Jan. 2005.
%
%**************************************************************************
% Author:           E. A. P. Habets, M. R. P. Thomas and A. H. Moore
% Date:             27 July 2010
% Version: $Id: sphBasis.m 9226 2016-12-21 12:09:45Z amoore1 $
% History:  2010-07-27: Original version
%           2016-12-21: Rewritten for speed - precompute factorials and
%                       scaling
%**************************************************************************

narginchk(3,3);

%enforce column vectors for az and inc
az = az(:);
inc = inc(:);

%some constants
nSH = (sphOrdMax+1)^2;                    %[1 1]; - total number of spherical harmonics
k_sphOrd = (2*(0:sphOrdMax).' +1)/(4*pi); %[sphOrdMax+1 1]; - use index sphOrd+1 to get sphOrd component
fact = factorial((0:2*sphOrdMax).');      %[2*sphOrdMax 1]; - use index i+1 to get factorial(i)

Y = zeros(length(az),nSH);  % preallocate matrix
Y(:,1) = sqrt(k_sphOrd(1)); % zeroth component is constant (omnidirectional)

for sphOrd = 1:sphOrdMax
    % N.B. legendre uses a recursive algorithm, so there is potential for
    % optimisation if all terms were returned from one call
    
    P = legendre(sphOrd,cos(inc)).';      %[nDirections, sphOrd+1]; - only defined for positive SH orders
                                          % - use index sphDeg+1 to get sphOrd component
    
    % Calculate 0 degree as special case
    col_i = sphOrd^2 + sphOrd + 1;        %use 0 degree index into output matrix as reference
    Y(:,col_i) = sqrt(k_sphOrd(sphOrd+1)) * P(:,1);
    
    % Loop over +ve degrees and obtain corresponding -ve degrees 
    for sphDeg = 1:sphOrd
        % C: normalisation factor
        % - indexes into the main precomputed normalisations k_sphOrd
        % - indexes into precomputed factorial functions
        C = sqrt(k_sphOrd(sphOrd+1)*fact(sphOrd-sphDeg+1)/fact(sphOrd+sphDeg+1)) ...
            * P(:,sphDeg+1);
        Y(:,col_i+sphDeg) = C .* exp(1i*sphDeg*az); %positive sphDeg;
        Y(:,col_i-sphDeg) = (-1)^-sphDeg .* conj(Y(:,col_i+sphDeg)); %negative sphDeg
    end
end
