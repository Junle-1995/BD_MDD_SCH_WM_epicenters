function [CenCoor] = gretna_centroid_coor_rois(Path_filename)

%==========================================================================
% This function is used to calculate centroid coordination of each ROI 
% included in the nifti image. The image is typically a parcellation atlas.
%
%
% Syntax: function [CenCoor] = gretna_centroid_coor_rois(Path_filename)
%
% Input:
%  Path_filename:
%                The directory & filename of the nifti image.
%
% Output:
%        CenCoor:
%                The centroid coordination of each region.
%
% Jinhui WANG, NKLCNL, BNU, BeiJing, 2014/09/13, jinhui.wang.1982@gmail.com
%==========================================================================

[pathstr, name, ext] = fileparts(Path_filename);

if ~isempty(pathstr)
    cd (pathstr)
end

Vout = spm_vol([name ext]);
[Ydata, ~] = spm_read_vols(Vout);

Nreg = max(Ydata(:));
voxcoor = ones(Nreg,4);

for i = 1:Nreg
    ind = find(Ydata == i);
    [I, J, K] = ind2sub(size(Ydata),ind);
    voxcoor(i,1) = mean(I);
    voxcoor(i,2) = mean(J);
    voxcoor(i,3) = mean(K);
end

CenCoor = Vout.mat*voxcoor';
CenCoor(4,:) = [];
CenCoor = CenCoor';

return