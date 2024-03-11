
%% BCI online
% run RDA based on exp type (train/test) and RDA type (unist/tsb) 

% 2021.09.03
% v1.0

%%

clear all; close all

mode = input('Play mode [train/test]: ','s');
rda = input('RDA type [unist/tsb]: ','s');

switch rda
    case 'unist'
        
        switch mode
            case 'train'
                RDA_train_unist();
            case 'test'
                RDA_test_unist();
                endsock
        end
        
    case 'tsb'
        
        switch mode
            case 'train'
                RDA_train_tsb();
            case 'test'
                RDA_test_tsb();
        end
end

                
 %% Get performance
cal_result('FAN'); % fan 
cal_result('AIRD'); % air dressor
cal_result('Def'); % lamp
cal_result('Def'); % doorlock
cal_result('Def'); % bluetooth speaker
cal_result('Def'); % air conditioner

    
