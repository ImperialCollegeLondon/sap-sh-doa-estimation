function[distmat] = az_inc_bin_distances(azg,incg)

ref_vec_inc = incg(:,1);
ref_vec_az = azg(:,1);
[ref_x,ref_y,ref_z] = mysph2cart(ref_vec_az,ref_vec_inc,ones(size(ref_vec_az)));
%off_vec_az = reshape(azg(:,1:180),[],1);
%off_vec_inc = reshape(incg(:,1:180),[],1);
off_vec_az = reshape(azg,[],1);
off_vec_inc = reshape(incg,[],1);
[off_x,off_y,off_z] = mysph2cart(off_vec_az,off_vec_inc,ones(size(off_vec_az)));
distmat = distcos([ref_x,ref_y,ref_z],[off_x,off_y,off_z],'x');