function EP = EpochData(sig,trig,param)

try
    switch param.EpocType
        case 'meantr'
            EP                      = Epoching(sig,trig,param);
            
        case 'meantr_tsb'
            EP                      = Epoching_TSBAR(sig, trig, param);
            
        case 'singletr'
            EP                      = Epoching_singletr(sig, trig, param);
            param.Target = EP.target;
            
        case 'singletr_tsb'
            EP                      = Epoching_singletr_TSBAR(sig, trig, param);
            param.Target = EP.target;
            
    end
catch
    keyboard;
end