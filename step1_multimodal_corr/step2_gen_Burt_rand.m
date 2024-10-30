%% run the Burt permutation
path = pwd;

py.importlib.import_module('Burt_rand')
load('Burt_para_fit.mat')

coord_file = [path '\Burt\cencoor.txt'];
data_path = [path '\Burt\data'];

Group = {'BD','MDD','SCH'};
Test_type = {'volume'};

for igroup = 1:length(Group)
    for itest = 1:length(Test_type)
        output_dir = [path '\Burt\data_rand\' Test_type{itest} '_' Group{igroup}];
        mkdir(output_dir)
        cd(output_dir)

        c = py.Burt_rand.volume_gen(coord_file,output_dir,...
                    [data_path '\' Test_type{itest} '_' Group{igroup} '.txt'],...
                    single(fit_best.([Test_type{itest} '_' Group{igroup}])(1)), ...
                    single(fit_best.([Test_type{itest} '_' Group{igroup}])(2)), ...
                    [Test_type{itest} '_' Group{igroup} '_burt_rand.mat']);
    end
end

cd(path)