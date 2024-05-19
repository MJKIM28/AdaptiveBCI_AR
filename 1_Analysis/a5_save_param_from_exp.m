clear all; close all;
%% save parameters

fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07',...
    'Subtest08','Subtest09','Subtest10','Subtest12','Subtest13','Subtest14','Subtest15'};

Nsub = length(SubNameList);


LAMBDA = 0.02;
%%
for s = 7:Nsub %subtest08 ~
    SubName = SubNameList{s};

%% parameters

    p = load([fpath,'\Dat_',SubName,'\param.mat']);
    UPmodel = p.param.update;

     save([fpath,'\Simulate_',num2str(LAMBDA),'param\',SubName],'UPmodel');
end
