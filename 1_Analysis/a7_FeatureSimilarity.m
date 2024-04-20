%% feature similarity 


clear all

LAMBDA = 0.02;
%%
fpath = 'E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot';
codepath = 'E:\1_연구\분석코드\2024_adaptive_BCI\AdaptiveBCI_AR\0_online_code';


SubNameList = {'Subtest02','Subtest03','Subtest04','Subtest05','Subtest07'};

Nsub = length(SubNameList);

Ntr_tr = 8;
Ntr_con = 15;
Nsess = 6;
Ntr_te = Ntr_con*Nsess;

Nch = 31;
Fs = 500;
Times = -0.2:1/Fs:0.6-1/Fs;
chanloc = readlocs('D:\[1]EEGBCI\[2]Research\Code&Algorithm\BP_channelLocs.locs');
%%
param.DSP.type = 'nf';
param.trD.mode = 'testing';
param.repeat = 10;
param.Stims = 1:4;
param.NumStims = length(param.Stims);
param.Fs = Fs;
param.Baseline = 0.2*param.Fs;
param.Epocline = 0.6*param.Fs;
param.Totalepoc = param.Baseline + param.Epocline;
%%
cd(codepath)
CS_all = []; CSused = []; NtrUsed = [];
MD_all = []; MDused = [];
for s = 1:Nsub
    SubName = SubNameList{s};
    load([fpath,'\Simulate_',num2str(LAMBDA),'param\',SubName])

    load([fpath,'\Test\Epoch\',SubName])

    load([fpath,'\Test\Block\',SubName])
    badch = Block.badch;
    chlist = setdiff(1:Nch,badch);
    chanloc_in = chanloc;
    chanloc_in(badch) = [];
    %%
    upcheck = [UPmodel.upcheck];
    uptrials = {UPmodel.upids};
    prediction = [UPmodel.prediction];
    target = [UPmodel.target];
    Correctness = prediction == target;
    mdl_init = UPmodel(1).mdl_init;
    DSP_init = UPmodel(1).DSP_init;


    FeatTar = []; FeatNar = [];
    for t = 1:length(Epoch) 


    
        if t <= 16
            DSPweights = DSP_init;
            MDL = mdl_init;
        else
            tr = t - 15;
             indices = find(upcheck(1:tr-1)==1);
             if ~isempty(indices)
             ind_use = indices(end);           
             DSPweights = UPmodel(ind_use).DSP;
             MDL = UPmodel(ind_use).mdl.mdl;
             else
                DSPweights = DSP_init;
                MDL = mdl_init;
             end

        end
        param.DSP = DSPweights;
        param.trD.mdl = MDL;
        ep = Epoch_condition(Epoch{t},param);

        [feat] = FeatureExt_DSP(ep,param);
        feat_re  = reshape(feat,param.repeat,param.NumStims,DSPweights.nf,[]);

        feat_tar = feat_re(:,ep.target,:,:);
        feat_nar = feat_re(:,setdiff(1:param.NumStims,ep.target),:,:);
     
        FeatTar(:,t,:,:) = feat_tar;
        FeatNar(:,:,t,:,:) = feat_nar;

    end

    feat_tar_pre = FeatTar(:,1:Ntr_con,:,:);
    feat_tar_main = FeatTar(:,Ntr_con+1:end-Ntr_con,:,:);
    feat_tar_post = FeatTar(:,end-Ntr_con+1:end,:,:);

    feattarpre_mean = squeeze(mean(reshape(feat_tar_pre,[],size(feat_tar_pre,3),size(feat_tar_pre,4)),'omitnan'))';
    feattarpre_mean = feattarpre_mean(:);


    feattarmain_temp = reshape(feat_tar_main,[],size(feat_tar_main,3),size(feat_tar_main,4));
    feattarmain_temp = permute(feattarmain_temp,[1,3,2]); % sample x time x comp
    feattarmain = reshape(feattarmain_temp,size(feattarmain_temp,1),[]);

    feattarpost_temp = reshape(feat_tar_post,[],size(feat_tar_main,3),size(feat_tar_main,4));
    feattarpost_temp = permute(feattarpost_temp,[1,3,2]); % sample x time x comp
    feattarpost = reshape(feattarpost_temp,size(feattarpost_temp,1),[]);

    feattarall = [feattarpre_mean';feattarmain;feattarpost];

    %-- cosine similarity
    cossim_temp = 1-pdist(feattarall,'cosine');
    cossim = squareform(cossim_temp);
    cossim_w_pre = cossim(1,2:end);
    CosSim_w_pre = reshape(cossim_w_pre,param.repeat,Ntr_con,[]);

    %-- all trial
    figure; imagesc(squeeze(mean(CosSim_w_pre)))

    %-- only used trials
    CStemp = reshape(CosSim_w_pre,param.repeat,[]);
    CSused = []; NtrUsed= [];
    for tt = 1:length(upcheck)
        if upcheck(tt)
        
            CSused=[CSused; CStemp(uptrials{tt},tt)];
            NtrUsed(tt) = length(uptrials{tt});
        end
    end
         CS_all(:,:,1:size(CosSim_w_pre,3),s) = CosSim_w_pre;
         CSused_all{s} = CSused;
         NtrUsed_all(1:size(CosSim_w_pre,3)*Ntr_con,s) = NtrUsed;


         %-- Mahalanobis distance

             feattarpre_tmp = reshape(feat_tar_pre,[],size(feat_tar_pre,3),size(feat_tar_pre,4));
    feattarpre_tmp = permute(feattarpre_tmp,[1,3,2]); % sample x time x comp
    feattarpre = reshape(feattarpre_tmp,size(feattarpre_tmp,1),[]);
    feattarpre_2 = feattarpre(:,1:60*2);

        feattarmain_2 = feattarmain(:,1:60*2);
        feattarpost_2 = feattarpost(:,1:60*2);

       MD =  sqrt(mahal(feattarmain_2,feattarpre_2));
       MD2 =  sqrt(mahal(feattarpost_2,feattarpre_2));
    MD_w_pre = reshape(MD,param.repeat,Ntr_con,[]);
    MD_w_pre2 = reshape(MD2,param.repeat,Ntr_con,[]);
    MD_w_pre_all = cat(3,MD_w_pre,MD_w_pre2);

    figure; imagesc(squeeze(mean(MD_w_pre_all)));
   
    %-- only used trials
    MDtemp = reshape(MD_w_pre_all,param.repeat,[]);
    MDused = []; 
    for tt = 1:length(upcheck)
        if upcheck(tt)
        
            MDused=[MDused; MDtemp(uptrials{tt},tt)];

        end
    end

 MD_all(:,:,1:size(MD_w_pre_all,3),s) = MD_w_pre_all;
         MDused_all{s} = MDused;



end