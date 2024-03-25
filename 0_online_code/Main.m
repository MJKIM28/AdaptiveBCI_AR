%% Experiment: Adaptive BCI in AR env
% 2024.03

% pre 1 session -> main 4 session -> post 2 session 

% Minju Kim


clear all; close all; clc

SubName = input('Subject (SubXX):','s');

% 15 trials x 7 sessions (1 pre, 4 main, 2 post)
Ntr = 15;
Nsess = 7;

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
fprintf('Ready for main exp\n')

pause;


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

intervals = [10 5 7 15 12 14 5 10 9 11 13 5 12 10; 9 13 5 8 7 10 12 6 10 15 9 12 6 10; 15 5 9 13 10 14 6 10 12 11 7 10 13 7];

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

MW = tcpclient('127.0.0.1',1800);
fopen(MW);


h = Figs.h;
axAll = Figs.axAll;
ax = Figs.ax;
stim_img = Figs.stim_img;


%% Pre (Control) 
fprintf('Ready for Pre session (without video)\n')
pause;

ClickCommand();
instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',3,54)
pause;

sess = 1;
set(h,'windowstate','minimized');
for tr = 1:Ntr
    fprintf('\nSession%d Trial %d\n..',sess,tr)
    fprintf('Stim start in %ds\n',3)
    pause(3)
    target = targetlist(tr,sess);
    targetpresent(h,MW,axAll,stim_img,target); % 3 sec

    %-- start BCI
    STARTBCI
    set(h,'windowstate','minimized');
    fwrite(MW,'57'); % trigger: start BCI
    fprintf('> BCI start\n')
    pause(15);
    fprintf('\nTrial %d end\n',tr)
end
fprintf('Pre session end\n')

%% Exp

fprintf('Run RDA_test_adaptive.exe and start stim with video\n')
pause;
for sess = 2:Nsess-2
    
    ClickCommand();
    instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',3,54)
    pause;
    set(h,'windowstate','minimized');

    for tr = 1:Ntr
        fprintf('\nSession%d Trial %d\n..',sess,tr)
        fprintf('Stim start in %ds\n',intervals(sess-1,tr))
        pause(intervals(sess-1,tr))
        target = targetlist(tr,sess);
        targetpresent(h,MW,axAll,stim_img,target); % 3 sec

        %-- start BCI
        STARTBCI
        set(h,'windowstate','minimized');
        fwrite(MW,'57'); % trigger: start BCI
        fprintf('> BCI start\n')
        pause(15);
        fprintf('\nTrial %d end\n',tr)
    end
end

%% post (control)
fprintf('Ready for Post session (without video)\n')
pause;

for sess = Nsess-1:Nsess
    
    ClickCommand();

    instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',3,54)
    pause;
    set(h,'windowstate','minimized');
    for tr = 1:Ntr
        fprintf('\nSession%d Trial %d\n..',sess,tr)
        fprintf('Stim start in %ds\n',3)
        pause(3)
        target = targetlist(tr,sess);
        targetpresent(h,MW,axAll,stim_img,target); % 3 sec

        %-- start BCI
        STARTBCI
        set(h,'windowstate','minimized');
        fwrite(MW,'57'); % trigger: start BCI
        fprintf('> BCI start\n')
        pause(15);
        fprintf('\nTrial %d end\n',tr)
    end
end
fprintf('Post session end\n')


%%
vars.target = targetlist;
vars.intervals = intervals;

save(SubName,'vars')

    instruction(Figs.h,MW,Figs.ax,'실험이 모두 종료되었습니다',10,54)
