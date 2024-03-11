function [EP] = Epoch_condition(EP, param)
Nt = size(EP.lat,3);

if strcmp(param.trD.mode,'training')
    
    STtar = [];
    STntar = [];
    
    for b = 1:Nt
        
        tr_tar = EP.lat(param.Baseline+1,:,b) == EP.target(b);
        tr_ntar = EP.lat(param.Baseline+1,:,b) ~= EP.target(b) & ismember(EP.lat(param.Baseline+1,:,b),param.Stims);
        tr_ntar(isnan(EP.lat(param.Baseline+1,:,b))) = 0 ; % ignore NaN
        
        STtar = cat(3,STtar,EP.dat(:,:,tr_tar,b));
        STntar = cat(3,STntar,EP.dat(:,:,tr_ntar,b));
        
    end
    basemean = mean(STtar(1:param.Baseline,:,:),1);
    STtar = STtar - repmat(basemean,param.Totalepoc,1,1);
    basemean = mean(STntar(1:param.Baseline,:,:),1);
    STntar = STntar - repmat(basemean,param.Totalepoc,1,1);
    
    EP.tar = STtar; % time x channel x (trial x block)
    EP.nar = STntar; % time x channel x (trial x block)
    
else
    ST = zeros(size(EP.dat,1),size(EP.dat,2),param.repeat,param.NumStims,Nt);
    for b = 1:Nt
        for n = 1:param.NumStims
            tr_tar = find(EP.lat(param.Baseline+1,:,b) == param.Stims(n));
            Nrepeat = length(tr_tar);
            if param.repeat ~= Nrepeat
                Nrepeat_new = min([Nrepeat param.repeat]);
                tr_tar = tr_tar(1:Nrepeat_new);
            else
                Nrepeat_new = Nrepeat;
            end
            ST(:,:,1:Nrepeat_new,n,b) = EP.dat(:,:,tr_tar,b);
            %     ST = cat(3,ST,EP.dat(:,:,tr_tar,b));
        end
    end
    basemean = mean(ST(1:param.Baseline,:,:,:,:),1);
    ST = ST - repmat(basemean,param.Totalepoc,1,1,1);
    
    EP.tar = [];
    EP.nar = ST; % time x channel x trial x stimulus x block
end