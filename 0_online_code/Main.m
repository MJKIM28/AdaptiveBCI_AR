%% Experiment: Adaptive BCI in AR env
% 2024.03

% pre 1 session -> main 4 session -> post 2 session 

% Minju Kim


clear all; close all; clc

SubName = input('Subject (SubXX):','s');

% 15 trials x 6 or 7 sessions (1 pre, 5 main, 1 or 2 post)
Ntr = 15;
Npost = 1;
Nsess = Npost + 6;


cd('C:\Users\minju\Desktop\mjkim\adaptive_BCI\AdaptiveBCI_AR\0_online_code')
%% Training
fprintf('Start Recording\n');
pause;    

fprintf('\nRun middleware (QtHue_v1.exe)\n');
fprintf('Hue address:192.168.1.231\n')
fprintf('Port: 4FF8\n')
pause;


fprintf('\nRun stimulus (Final1.exe)\n')
fprintf('Select random sequence\n')
fprintf('Select Light\n')
fprintf('Connect to server\n')
pause;

run('RDA_train_adaptive.m')

%% Ready for Exp


screens = get(0,'MonitorPositions');
fprintf('%d screen detected\n',size(screens,1));
ScreenSelect = input('Which screen?');
screen_use = screens(ScreenSelect,:);

h = figure(12);
set(h,'position',screen_use ...
    ,'WindowState','fullscreen','menubar','none','toolbar','none','color',[0 0 0])
winsize = get(h,'Position');
% winsize = [1 1 1920 1080];

set(h,'WindowState','minimized')
Figs.h = h;

intervals = [10 5 7 15 12 14 5 10 9 11 13 5 12 10 11; ...
    9 13 5 8 7 10 12 6 10 15 9 12 6 10 12; ...
    15 5 9 13 10 14 6 10 12 11 7 10 13 7 10;...
    7 8 5 10 15 9 13 12 8 14 7 9 10 11 6;
    9 13 10 14 8 15 12 7 11 14 13 10 14 15 6];

stims = [1 2 3 4];
NSeq = Ntr*Nsess;
targets = repmat(stims, 1,ceil(NSeq/length(stims)));

targets_ind = randperm(length(targets));
targetsTemp = targets(targets_ind(1:NSeq));

target_all = targetsTemp(:);
same_id = find(diff(target_all)==0) +1;
for ii = 1:length(same_id)

    while target_all(same_id(ii)) ~= 4 ... % 4-> color change
            && (target_all(same_id(ii)) == target_all(same_id(ii)-1))
        if same_id(ii) < length(target_all)
            potential = setdiff(stims,[target_all(same_id(ii)-1) target_all(same_id(ii)+1)]);
        else
            potential = setdiff(stims,target_all(same_id(ii)-1));
        end
        target_all(same_id(ii)) =potential(randi(length(potential)));

    end
end
targetlist = reshape(target_all,Ntr,Nsess);


stim_img = cell(4,2);
stim_img{1,1} = imread('./image/lightfullon.png');
stim_img{2,1} = imread('./image/lighthalfon.png');
stim_img{3,1} = imread('./image/lightonoff.png');
stim_img{4,1} = imread('./image/color.png');
stim_img{1,2} = imread('./image/lightfullon_target.png');
stim_img{2,2} = imread('./image/lighthalfon_target.png');
stim_img{3,2} = imread('./image/lightonoff_target.png');
stim_img{4,2} = imread('./image/color_target.png');
Figs.stim_img = stim_img;


Figs.axAll = targetmake(Figs.h,winsize);
Figs.ax = makeinstruction(Figs.h,winsize);


%% Pre (Control) 
fprintf('Ready for Pre session (without video)\n')
fprintf('Run RDA_test_adaptive.exe and start stim without video\n')

pause;

MW = tcpclient('127.0.0.1',1800);
fopen(MW);


instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',3,54)
ClickCommand(1300,950);
pause;

sess = 1;
set(Figs.h,'windowstate','minimized');
for tr = 1:Ntr
    fprintf('\nSession%d Trial %d\n..',sess,tr)
    fprintf('Stim start in %ds\n',3)
    pause(3)
    target = targetlist(tr,sess);
    targetpresent(Figs.h,MW,Figs.axAll,Figs.stim_img,target); % 3 sec

    %-- start BCI
    STARTBCI
    set(Figs.h,'windowstate','minimized');
    fwrite(MW,'57'); % trigger: start BCI
    fprintf('> BCI start\n')
    pause(15);
    fprintf('\nTrial %d end\n',tr)
end
fprintf('Pre session end\n')
clear MW

try
    pause;
load(['Dat_',SubName,'/param.mat']);

pred_pre =  param.prediction(end-Ntr+1:end);
vars.Acc_pre = mean(pred_pre' == targetlist(:,1));
fprintf('==============\nAcc_pre  %.2f \n==============\n',vars.Acc_pre);
catch
    fprintf('!! Error in accuracy computation\n')
end

%% Exp
fprintf('Ready for Main session (with video)\n')
fprintf('Run RDA_test_adaptive.exe and start stim with video\n')
pause;

MW = tcpclient('127.0.0.1',1800);
fopen(MW);

for sess = 2:Nsess-Npost
    
    ClickCommand(1300,950);
    instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',3,54)
    pause;
    set(Figs.h,'windowstate','minimized');

    for tr = 1:Ntr
        fprintf('\nSession%d Trial %d\n..',sess,tr)
        fprintf('Stim start in %ds\n',intervals(sess-1,tr))
        pause(intervals(sess-1,tr))
        target = targetlist(tr,sess);
        targetpresent(Figs.h,MW,Figs.axAll,Figs.stim_img,target); % 3 sec

        %-- start BCI
        STARTBCI
        set(Figs.h,'windowstate','minimized');
        fwrite(MW,'57'); % trigger: start BCI
        fprintf('> BCI start\n')
        pause(18);
        fprintf('\nTrial %d end\n',tr)
    end
end
clear MW

try
    pause;
load(['Dat_',SubName,'/param.mat']);

pred_main = param.prediction(end-(Nsess-2)*Ntr+1:end);
target_main =  targetlist(:,2:6);
vars.Acc_mainAll = mean(pred_main'==target_main(:));

for s= 1:Nsess-1-Npost
    vars.Acc_main_sess(s) = mean(pred_main(Ntr*(s-1)+1:s*Ntr)' == targetlist(:,1+s));
end
fprintf('==============\nAcc_main1  %.2f \nAcc_main2  %.2f\nAcc_main3  %.2f\nAcc_main4  %.2f\nAcc_main(all) %.2f\n====================\n',...
    vars.Acc_main_sess(1),vars.Acc_main_sess(2),vars.Acc_main_sess(3),vars.Acc_main_sess(4),vars.Acc_mainAll);
catch
    fprintf('!! Error in accuracy computation\n')
end

%% post (control)
fprintf('Ready for Post session (without video)\n')
fprintf('Run RDA_test_adaptive.exe and start stim without video\n')
pause;

MW = tcpclient('127.0.0.1',1800);
fopen(MW);
for sess = Nsess-Npost+1:Nsess
    
    ClickCommand(1300,950);
    instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',3,54)
    pause;
    set(Figs.h,'windowstate','minimized');
    for tr = 1:Ntr
        fprintf('\nSession%d Trial %d\n..',sess,tr)
        fprintf('Stim start in %ds\n',5)
        pause(3)
        target = targetlist(tr,sess);
        targetpresent(Figs.h,MW,Figs.axAll,Figs.stim_img,target); % 3 sec

        %-- start BCI
        STARTBCI
        set(Figs.h,'windowstate','minimized');
        fwrite(MW,'57'); % trigger: start BCI
        fprintf('> BCI start\n')
        pause(18);
        fprintf('\nTrial %d end\n',tr)
    end
end
fprintf('Post session end\n')

try
load(['Dat_',SubName,'/param.mat']);

pred_post = param.prediction(end-(Nsess-6)*Ntr+1:end);
target_post = targetlist(:,7:end);
vars.Acc_post = mean(pred_post'==target_post(:));
fprintf('==============\nAcc_post %.2f\n',...
    vars.Acc_post);

catch
    fprintf('!! Error in accuracy computation\n')
end
%%
instruction(Figs.h,MW,Figs.ax,'실험이 모두 종료되었습니다',10,54)


vars.target = targetlist;
vars.intervals = intervals;
vars.subname = SubName;

%-- Accuracy
load(['Dat_',SubName,'/param.mat']);
close
pred_pre =  param.prediction(end-Nsess*Ntr+1:end-(Nsess-1)*Ntr);
pred_main = param.prediction(end-(Nsess-1)*Ntr+1:end-(Nsess-6)*Ntr);
pred_post = param.prediction(end-(Nsess-6)*Ntr+1:end);

vars.Acc_pre = mean(pred_pre' == targetlist(:,1));

target_main =  targetlist(:,2:6);
vars.Acc_mainAll = mean(pred_main'==target_main(:));

target_post = targetlist(:,7:end);
vars.Acc_post = mean(pred_post'==target_post(:));

for s= 1:Nsess-1-Npost
    vars.Acc_main_sess(s) = mean(pred_main(Ntr*(s-1)+1:s*Ntr)' == targetlist(:,1+s));
end
fprintf('==============\nAcc_pre  %.2f \n \nAcc_main1  %.2f \nAcc_main2  %.2f\nAcc_main3  %.2f\nAcc_main4  %.2f\nAcc_main(all) %.2f\n\nAcc_post %.2f\n',...
    vars.Acc_pre,vars.Acc_main_sess(1),vars.Acc_main_sess(2),vars.Acc_main_sess(3),vars.Acc_main_sess(4),vars.Acc_mainAll,vars.Acc_post);

%-- save variables
save(SubName,'vars')
