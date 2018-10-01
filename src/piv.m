function[I,extras] = piv(P)
%
%
% Inputs:
%       P compensated eigenbeams of order 0 an 1 [P00, P1(-1) P10 P11]


P = P(:,1:4,:); % discard uneeded eigenbeams

% monopole is just 0-order eigenbeam
% dipole is created by steering PWD beams towards negative x,y,z directions
steering_angles     = [ pi,    pi/2 ;...             % [azimuth, inclination]
                       -pi/2,  pi/2;...
                        0,     pi];
ndipoles = size(steering_angles,1);                   
B = sphBasis(steering_angles(:,1),steering_angles(:,2),1);                 % [ndipoles, n_eigenbeams]
B(:,1) = [];
ord_1_eigenbeams = P(:,2:4,:);
[nfreqs,neigbeams,nframes] = size(ord_1_eigenbeams);
dipoles = reshape(permute(ord_1_eigenbeams,[1 3 2]),nfreqs*nframes,neigbeams) * B.';	   % [nfreqs*nframes, ndipoles]
dipoles = permute(reshape(dipoles,nfreqs,nframes,ndipoles),[1 3 2]);       % [nfreqs,ndipoles,nframes]


% intensity
I = -0.5*real( bsxfun(@times,conj(P(:,1,:)),dipoles));                     % [nfreqs, naxes, nframes] (naxes: x,y,z)

if nargout>1
    extras.monopole = P(:,1,:);
    extras.dipoles = dipoles;
end