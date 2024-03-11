function [EP_new, param] = FeatExt_basic_prep(EP, param)

Nb = size(EP.dat,4);

EP_stim = NaN(param.Totalepoc,param.NumCh,param.NumStims,Nb);

for s = 1:param.NumStims
    EP_stim_lat = EP.lat(param.Baseline+1,:,:) == s;
    for b = 1:Nb
        % average repetitions
        EP_stim(:,:,s,b) = mean(EP.dat(:,:,EP_stim_lat(1,:,b),b),3);
    end
end

% baseline correction
basemean = mean(EP_stim(1:param.Baseline,:,:,:),1);
EP_stim = EP_stim - repmat(basemean, param.Totalepoc,1,1,1);

EP_new.tar                 = zeros(param.NumCh, param.Totalepoc, Nb);
if strcmp(param.trD.mode,'training')
    EP_new.nar                 = zeros(param.NumCh, param.Totalepoc, param.NumStims-1, Nb);
else
    EP_new.nar                 = zeros(param.NumCh, param.Totalepoc, param.NumStims, Nb);
end

for b = 1:Nb
    if strcmp(param.trD.mode,'training')
        % separate target and nontarget trials for every single block
        EP_new.tar(:,:,b) = permute(EP_stim(:,:,EP.target(b),b),[2,1,3]);
        EP_new.nar(:,:,:,b) = permute(EP_stim(:,:,setdiff(param.Stims,EP.target(b)),b),[2,1,3,4]);
    else
        EP_new.nar(:,:,:,b) = permute(EP_stim,[2,1,3,4]);
    end
end

EP_new.badch = param.badch;
EP_new.target = EP.target;

for k = 1:param.NumCh
    set(param.h(k,1), 'XData', param.Time,'YData', squeeze(mean(EP_new.tar(k,:,:),3)), 'Color','r','LineWidth',2);
    set(param.h(k,2), 'XData', param.Time,'YData', squeeze(mean(EP_new.nar(k,:,1,:),4)),'Color','g');
    set(param.h(k,3), 'XData', param.Time,'YData', squeeze(mean(EP_new.nar(k,:,2,:),4)),'Color','b');
    set(param.h(k,4), 'XData', param.Time,'YData', squeeze(mean(EP_new.nar(k,:,3,:),4)),'Color','k');
    ylim(param.SH(k),[-2 2]);
    if ~strcmp(param.decoder.mode,'training')
        set(param.h(k,1), 'XData', param.Time,'YData', squeeze(mean(EP_new.nar(k,:,4,:),4)),'Color','m');
           ylim(param.SH(k),[-5 5]);
           legend({'','1','2','3','4'});
    end
   title(param.SH(k),param.Ch{k});
   xlim(param.SH(k),[-0.2 0.6])
end



