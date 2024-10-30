%% try all the para to Burt permutation
%--------------------------------------------------------------------------
% install python3.7/3.8
% download and install brainsmash (https://github.com/murraylab/brainsmash)
%
% edit the file (mapgen\eval.py):
%   add: "from scipy import stats" at the first line
%   add: "import matplotlib" and "matplotlib.use('Agg')" before: "import matplotlib.pyplot as plt"
%   add: "rho = scipy.stats.spearmanr(emp_var,mu)[0]" and "return fig, rho" at the last line
%
%--------------------------------------------------------------------------
%% para
path = pwd;

Group = {'BD','MDD','SCH'};
Test_type = {'volume'};


%% real data
data_path = [path '\Burt\data'];
mkdir(data_path);

% signal
cd(path)
cd ..
load('Results_signal.mat')
for igroup = 1:length(Group)
    Z_file = [path '\Burt\data\volume_' Group{igroup} '.txt'];
    fid = fopen(Z_file,'wt');
    for i = 1:size(results.volume_Z,1)
        fprintf(fid, '%f\n', results.volume_Z(i,igroup));
    end
    fclose(fid);
end


%% get coordinate
cd(path)
cd ..
[CenCoor] = gretna_centroid_coor_rois([pwd '\atlas\JHU-WhiteMatter-labels_1.5mm.nii']);
cd(path)

coord_file = [path '\Burt\cencoor.txt'];
fid = fopen(coord_file,'wt');
for i = 1:size(CenCoor,1)
    fprintf(fid, '%f %f %f\n', CenCoor(i,:));
end
fclose(fid);


%% try para
py.importlib.import_module('Burt_rand_para_set')

fit_results = struct;
fit_best = struct;
ns = 16:4:48;
knn = 16:4:48;

for igroup = 1:length(Group)
    for itest = 1:length(Test_type)
        output_dir = [path '\Burt\para_pic\' Test_type{itest} '_' Group{igroup}];
        mkdir(output_dir)
        cd(output_dir)

        disp(['Now processing the data in ' Test_type{itest} '_' Group{igroup} '  |' datestr(clock)])
        fit_results.([Test_type{itest} '_' Group{igroup}]) = [];
        for i = ns
            for j = knn
                pic_name = [num2str(i) '_' num2str(j)];

                c = py.Burt_rand_para_set.volume_gen(coord_file,output_dir,...
                    [data_path '\' Test_type{itest} '_' Group{igroup} '.txt'],...
                    single(i),single(j),pic_name);

                pic_file = dir([num2str(i) '_' num2str(j) '*']);
                if length(pic_file) ~= 1, error('?'), end
                rho = pic_file(1).name;
                rho_posi = strfind(rho,'_');
                rho = str2double(rho(rho_posi(2)+1:end-4));

                fit_results.([Test_type{itest} '_' Group{igroup}])(end+1,:) = [i,j,rho];
            end
        end
        [~,I] = max(fit_results.([Test_type{itest} '_' Group{igroup}])(:,3));
        fit_best.([Test_type{itest} '_' Group{igroup}]) = fit_results.([Test_type{itest} '_' Group{igroup}])(I,:);
    end
end

cd(path)
save('Burt_para_fit.mat','fit_results','fit_best')