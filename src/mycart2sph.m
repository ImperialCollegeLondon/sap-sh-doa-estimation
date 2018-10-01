function [az,inc,r] = mycart2sph(x,y,z)

%   3D Cartesian to Spherical Coordinates
%
%   [az,inc,r] = mycart2sph(x,y,z)
%
%   Inputs:
%       x       x-coordinate
%       y       y-coordinate
%       z       z-coordinate
%
%   Outputs:
%       az      Azimuth angle in radians, counterclockwise from xy plane 
%               from the positive x axis (otherwise referred to as phi)
%       inc     Inclination angle in radians, from positive z axis  
%				(otherwise referred to as theta) 
%       r       Radius
%
%   Notes:
%       The MATLAB function cart2sph reverses phi and theta.
%
%**************************************************************************
% Author:           E. A. P. Habets and M. R. P. Thomas and D. P. Jarrett
% Date:             24 August 2011
%**************************************************************************

if nargin==1
    [~,j,k] = size(x);
    if k==1 && j== 3   
        z = x(:,3);
        y = x(:,2);
        x = x(:,1);
    else
        error('Inputs must be 3 column vectors or a Mx3 matrix')
    end
end

hypotxy = hypot(x,y);
r = hypot(hypotxy,z);
inc = acos(z./r);
az = atan2(y,x);
