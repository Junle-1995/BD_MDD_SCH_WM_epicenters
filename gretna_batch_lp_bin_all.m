function [Lp_group, Lp_individual] = gretna_batch_lp_bin_all(Mats, Thr1, Thr2, Delta, Stype, Ttype)

%==========================================================================
% This function is used to calculate several global and nodal network
% metrics for both binary and weighted networks over a continuous threshold
% range for multiple subjects.
%
%
% Syntax: function [Net_para, Node_para] = gretna_batch_networkanalysis(Mat_path, File_filter, Thr1, Thr2, Delta, N_rand, Stype, Ttype)
%
% Inputs:
%          Mat_path:
%                   The directory where individual matrices are sorted.
%       File_filter:
%                   The prefix of those matrices.
%              Thr1:
%                   The lower bound of the threshold range.
%              Thr2:
%                   The upper bound of the threhsold range.
%             Delta:
%                   The interval of the threshold range.
%            N_rand:
%                   The number of random networks.
%             Stype:
%                   'pos' for only positive elements;
%                   'neg' for only negative elements using absolute values;
%                   'abs' for all elements using absolute values.
%             Ttype:
%                   'r' for correlation threshold;
%                   's' for sparsity threshold;
%                   'k' for edge threshold.
%
%
% Jinhui WANG, IBRR, SCNU, Guangzhou, 2019/10/31, jinhui.wang.1982@gmail.com
%==========================================================================

Thres = Thr1:Delta:Thr2;

Lp = zeros(length(Thres),size(Mats,3),size(Mats,1),size(Mats,1));
Lp_individual = zeros(size(Mats,3),size(Mats,1),size(Mats,1));

for isub = 1:size(Mats,3) % subjects
    disp(['Calculating network parameters for subject ' num2str(isub) '/' num2str(size(Mats,3)) '  |' datestr(clock)]);
    
    Matrix = Mats(:,:,isub);

    check_sym = Matrix - Matrix';
    if min(abs(check_sym(:))) > 10*eps
        error('?')
    else
        Matrix = (Matrix+Matrix')/2;
    end
    
    for ithres = 1:length(Thres) % thresholds
        [Bin, ~] = gretna_R2b_MST(Matrix,Stype,Ttype,Thres(ithres));

        D = gretna_distance(Bin);

        D(1:length(D)+1:end) = nan;
        
        Lp(ithres,isub,:,:) = D;
        
    end % thresholds

    for inode = 1:size(Mats,1)
        Lp_individual(isub,inode,:)   = gretna_auc(squeeze(Lp(:,isub,inode,:)),Delta);
    end
    
end % subjects

Lp_group = squeeze(mean(Lp_individual,1));

return