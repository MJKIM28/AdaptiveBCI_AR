%% ERP r2 (correlation with label)

clear all; close all
LAMBDA = 0.02;
%% parameters
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';

SubName = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07'};
Nsub = length(SubName);

Fs = 500;
Nch = 31;
Tbase = 0.2;
Tepoc = 0.6;
Lepoc = (Tbase+Tepoc)*Fs;

Ntr_tr = 8;
Ntr_con  = 15;
Nsess_main = 5;

Nrepeat = 10;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');

ChList = {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
Times = -0.2:1/Fs:0.6-1/Fs;

%% load (Epoch)
EP_tar_tr = NaN(Lepoc,Nch,Nrepeat*Ntr_tr,Nsub);
EP_tar_pre = NaN(Lepoc,Nch,Nrepeat*Ntr_con,Nsub);
EP_tar_main = NaN(Lepoc,Nch,Nrepeat*Ntr_con,Nsess_main,Nsub);
EP_tar_post = NaN(Lepoc,Nch,Nrepeat*Ntr_con,Nsub);
EP_nar_tr = NaN(Lepoc,Nch,Nrepeat*Ntr_tr*3,Nsub);
EP_nar_pre = NaN(Lepoc,Nch,Nrepeat*Ntr_con*3,Nsub);
EP_nar_main = NaN(Lepoc,Nch,Nrepeat*Ntr_con*3,Nsess_main,Nsub);
EP_nar_post = NaN(Lepoc,Nch,Nrepeat*Ntr_con*3,Nsub);
for s = 1:Nsub

    e_tr = load([fpath,'\Train\Epoch_15HZLPF\',SubName{s}]);
    e_te = load([fpath,'\Test\Epoch_15HZLPF\',SubName{s}]);
    chlist = setdiff(1:Nch,e_tr.Epoch.badch);

    Nsess_main_in = fix(length(e_te.Epoch) - Ntr_con*2)./Ntr_con;
    %-- target
    EP_tar_tr(:,chlist,:,s) = e_tr.Epoch.tar;
    EP_tmp1 = [];
    for tr = 1:Ntr_con
        EP_tmp1 = cat(3,EP_tmp1,e_te.Epoch{tr}.tar);
    end
    EP_tmp2 = [];
    for tr = Ntr_con+1:Ntr_con*(Nsess_main_in+1)
        EP_tmp2 = cat(3,EP_tmp2,e_te.Epoch{tr}.tar);
    end
    EP_tmp3 = [];
    for tr = Ntr_con*(Nsess_main_in+1)+1:Ntr_con*(Nsess_main_in+2)
        EP_tmp3 = cat(3,EP_tmp3,e_te.Epoch{tr}.tar);
    end
    EP_tar_pre(:,chlist,:,s) = EP_tmp1;
    EP_tar_main(:,chlist,:,1:Nsess_main_in,s) = reshape(EP_tmp2,[Lepoc,length(chlist),Nrepeat*Ntr_con,Nsess_main_in]);
    EP_tar_post(:,chlist,:,s) = EP_tmp3;


    EP_nar_tr(:,chlist,:,s) = e_tr.Epoch.nar;
    EP_nar_tmp1 = [];
    for tr = 1:Ntr_con
        EP_nar_tmp1 = cat(3,EP_nar_tmp1,e_te.Epoch{tr}.nar(:,:,1:3*Nrepeat));
    end
    EP_nar_tmp2 = [];
    for tr = Ntr_con+1:Ntr_con*(Nsess_main_in+1)
        EP_nar_tmp2 = cat(3,EP_nar_tmp2,e_te.Epoch{tr}.nar);
    end
    EP_nar_tmp3 = [];
    for tr = Ntr_con*(Nsess_main_in+1)+1:Ntr_con*(Nsess_main_in+2)
        EP_nar_tmp3 = cat(3,EP_nar_tmp3,e_te.Epoch{tr}.nar);
    end
    EP_nar_pre(:,chlist,:,s) = EP_nar_tmp1;
    EP_nar_main(:,chlist,:,1:Nsess_main_in,s) = reshape(EP_nar_tmp2,[Lepoc,length(chlist),Nrepeat*Ntr_con*3,Nsess_main_in]);
    EP_nar_post(:,chlist,:,s) = EP_nar_tmp3;


end
%% r2 - all trials

R2_pre = []; R2_main = []; R2_post = [];
for s= 1:Nsub

    r_pre = getR(EP_tar_pre(:,:,:,s),EP_nar_pre(:,:,:,s));
    r_main = [];
    for c = 1:size(EP_tar_main,4)
        r_main_temp = getR(EP_tar_main(:,:,:,c,s),EP_nar_main(:,:,:,c,s));
        if ~isempty(r_main_temp)
            r_main(:,:,c) = r_main_temp;
        end
    end
    r_post = getR(EP_tar_post(:,:,:,s),EP_nar_post(:,:,:,s));

    R2_pre(:,:,s) = r_pre.^2;
    R2_main(:,:,1:size(r_main,3),s) = r_main.^2;

    R2_post(:,:,s) = r_post.^2;

end


%% r2 - used trials
EP_tar_main_trial = reshape(EP_tar_main,size(EP_tar_main,1),size(EP_tar_main,2),Nrepeat,Ntr_con*size(EP_tar_main,4),Nsub);
EP_nar_main_trial = reshape(EP_nar_main,size(EP_tar_main,1),size(EP_tar_main,2),3,Nrepeat,Ntr_con*size(EP_tar_main,4),Nsub);
EP_tar_post_trial = reshape(EP_tar_post,size(EP_tar_post,1),size(EP_tar_post,2),Nrepeat,Ntr_con,Nsub);
EP_nar_post_trial = reshape(EP_nar_post,size(EP_tar_post,1),size(EP_tar_post,2),3,Nrepeat,Ntr_con,Nsub);

R2_main_partial = []; R2_post_partial = [];
for s= 1:Nsub
    load([fpath,'\Simulate_',num2str(LAMBDA),'param\',SubName{s}])
    upcheck = [UPmodel.upcheck];
    uptrial = {UPmodel.upids};

    ep_tar = []; ep_nar = [];
    for tr = 1:length(upcheck)-Ntr_con
        if upcheck(tr)
            ep_tar = cat(3,ep_tar, EP_tar_main_trial(:,:,uptrial{tr},tr,s));
            ep_nar = cat(4,ep_nar, EP_nar_main_trial(:,:,:,uptrial{tr},tr,s));
        end
    end
    ep_nar2 = reshape(ep_nar,size(ep_nar,1),size(ep_nar,2),[]);

    r_main = getR(ep_tar,ep_nar2);


    ep_tar = []; ep_nar = [];
    tt = 1;
    for tr = length(upcheck)-Ntr_con+1:length(upcheck)
        if upcheck(tr)
            ep_tar = cat(3,ep_tar, EP_tar_post_trial(:,:,uptrial{tr},tt,s));
            ep_nar = cat(4,ep_nar, EP_nar_post_trial(:,:,:,uptrial{tr},tt,s));
        end
        tt = tt+1;
    end
    ep_nar2 = reshape(ep_nar,size(ep_nar,1),size(ep_nar,2),[]);

    r_post = getR(ep_tar,ep_nar2);

    R2_main_partial(:,:,s) = r_main.^2;
    R2_post_partial(:,:,s) = r_post.^2;

end
%%
figure;
for s= 1:Nsub
    subplot(Nsub,7,(s-1)*7+1);
    imagesc(R2_pre(:,:,s)');
    %     caxis([0 0.08]);
    colormap jet
    colorbar

    for c = 1:size(R2_main,3)
        subplot(Nsub,7,(s-1)*7+1+c);
        imagesc(R2_main(:,:,c,s)');
        %     caxis([0 0.08]);
        colormap jet
        colorbar
    end

    subplot(Nsub,7,(s-1)*7+7);
    imagesc(R2_post(:,:,s)');
    %     caxis([0 0.08]);
    colormap jet
    colorbar


end

figure;
for s= 1:Nsub
    subplot(Nsub,2,(s-1)*2+1);
    imagesc(R2_main_partial(:,:,s)');
    %     caxis([0 0.08]);
    colormap jet
    colorbar

    subplot(Nsub,2,(s-1)*2+2);
    imagesc(R2_post_partial(:,:,s)');
    %     caxis([0 0.08]);
    colormap jet
    colorbar


end

%%
function r = getR(tar,nar)
% tar: time x ch x trial
% nar: time x ch x trial
EPs = cat(3,tar,nar);

labels = [ones(size(tar,3),1);zeros(size(nar,3),1)];

EPs_re = reshape(EPs,size(EPs,1)*size(EPs,2),[])';
r = [];
for tt = 1:size(EPs_re,2)
    xx = EPs_re(:,tt);
    yy = labels;
    rem = find(isnan(xx));
    xx(rem) = []; yy(rem) = [];
    if ~isempty(yy)
        r(tt) = pointbiserial(yy,xx);
    else
        r(tt) = NaN;
    end

end
r = reshape(r, size(EPs,1),size(EPs,2));


end

%% 3. amplitude (box plot)

%% 4. latency (box plot)

%% Statistic
%-- Friedman
%% 1. amplitude

%% 2. latency




