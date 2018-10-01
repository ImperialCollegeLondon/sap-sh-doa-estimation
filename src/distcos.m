function d=distcos(x,y,mode,w)
%DISTCOS calculate cosine distance D=(X,Y,MODE,W)
%
% Inputs: X,Y         Vector sets to be compared. Each row contains a data vector.
%                     X and Y must have the same number of columns.
%
%         MODE        Character string selecting the following options:
%                         'x'  Calculate the full distance matrix from every row of X to every row of Y
%                         'd'  Calculate only the distance between corresponding rows of X and Y
%                              The default is 'd' if X and Y have the same number of rows otherwise 'x'.
%
%         W           Optional weighting matrix: the distance calculated is (x-y)*W*(x-y)'
%                     If W is a vector, then the matrix diag(W) is used.
%
% Output: D           If MODE='d' then D is a column vector with the same number of rows as the shorter of X and Y.
%                     If MODE='x' then D is a matrix with the same number of rows as X and the same number of columns as Y'.
%

%      Copyright (C) Mike Brookes 1998, Hacked by Alastair Moore, April
%      2014
%      Version: $Id: disteusq.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx,p]=size(x); ny=size(y,1);
x = double(x); %AHM ensure double precision here (rather than on every usage)
y = double(y);

if nargin<3 | isempty(mode) mode='0'; end
if any(mode=='d') | (mode~='x' & nx==ny)

    % Do pairwise distance calculation

    nx=min(nx,ny);
    z=sum(x(1:nx,:) .* y(1:nx,:),2) ./ (cartnorm(x(1:nx,:)) .* cartnorm(x(1:nx,:)));
%     if nargin<4
%         d=sum(z.*conj(z),2);
%     else
%         error('''w'' parameter not supported')
    %elseif min(size(w))==1
%         
%         wv=w(:).';
%         d=sum(z.*wv(ones(size(z,1),1),:).*conj(z),2);
%     else
%         d=sum(z*w.*conj(z),2);
%     end
else
    
    % Calculate full distance matrix
    
    if p>1
        
        % x and y are matrices
        
        if nargin<4
            nxy=nx*ny;
            X = reshape(permute(x(:,:,ones(1,ny)),[1 3 2]),nxy,p); %repeat each value in x ny times in 3rd dimension then rearrange dimensions using permute such that a reshape ends up with 2nd dimension back in 2nd dimension 
            Y = reshape(permute(y(:,:,ones(1,nx)),[3 1 2]),nxy,p); %as above but the permute has 1and3 swapped so that we get all combinations
            %i.e. after rejigging 
            % X is [x(1,:);x(2,:);x(3,:);....;x(1,:);x(2,:);x(3,:);...]
            % y is [y(1,:);y(1,:);y(1,:);....;y(2,:);y(2,:);y(2,:);...]
            z = sum(X .* Y,2) ./ (cartnorm(X) .* cartnorm(Y));
            %now put vector into matrix so all (ix,iy) coordinates
            %represent the indices of original matrices which distances
            %relate to
            z = reshape(z,nx,ny);
            %z=permute(double(x(:,:,ones(1,ny))),[1 3 2])-permute(double(y(:,:,ones(1,nx))),[3 1 2]);
            %d=sum(z.*conj(z),3);
%         else
%             nxy=nx*ny;
%             z=reshape(permute(double(x(:,:,ones(1,ny))),[1 3 2])-permute(double(y(:,:,ones(1,nx))),[3 1 2]),nxy,p);
%             if min(size(w))==1
%                 wv=w(:).';
%                 d=reshape(sum(z.*wv(ones(nxy,1),:).*conj(z),2),nx,ny);
%             else
%                 d=reshape(sum(z*w.*conj(z),2),nx,ny);
%             end
%         end
%     else
%         
%         % x and y are vectors
%         
%         z=double(x(:,ones(1,ny)))-double(y(:,ones(1,nx))).';
%         if nargin<4
%             d=z.*conj(z);
%         else
%             d=w*z.*conj(z);
        end
    end
end
% if any(mode=='s')
%     d=sqrt(d);
% end
%d=abs(z);
%d = 1-acos(z)/pi;
d = acos(z)/pi;
