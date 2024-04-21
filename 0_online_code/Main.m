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


% cd('C:\Users\minju\Desktop\mjkim\adaptive_BCI\AdaptiveBCI_AR\0_online_code')
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

%-- set screen
h = figure(12);
set(h,'position',screen_use ...
    ,'WindowState','fullscreen','menubar','none','toolbar','none','color',[0 0 0])
winsize = get(h,'Position');
% winsize = [1 1 1920 1080];

set(h,'WindowState','minimized')
Figs.h = h;

h2 = figure(13);
set(h2,'position',screen_use ...
    ,'WindowState','fullscreen','menubar','none','toolbar','none','color',[0 0 0])
winsize = get(h2,'Position');
% winsize = [1 1 1920 1080];

set(h2,'WindowState','minimized')
Figs.h2 = h2;

h3 = figure(14);
set(h3,'position',screen_use ...
    ,'WindowState','fullscreen','menubar','none','toolbar','none','color',[0 0 0])
winsize = get(h3,'Position');
% winsize = [1 1 1920 1080];

set(h3,'WindowState','minimized')
Figs.h3 = h3;

Figs.axAll = targetmake(Figs.h,winsize);
Figs.ax = makeinstruction(Figs.h,winsize);
Figs.ax_h3 = makeinstruction(Figs.h3,winsize);


%-- NASA TLX
loadtypes = {'정신적 난이도 (Mental Demand)';'시간적 난이도 (Temporal Demand)';'신체적 난이도 (Physical Demand)';...
    '성취감 (Performance)';'노력 수준 (Effort)';'좌절감 (Frustration)'};
descriptions = {['작업 수행시 정신적/인지적 활동이 요구된 정도. 작업이 쉬웠나요, 어려웠나요? 단순했나요, 복잡했나요?',...
    newline,'How much mental and perceptual activity was required? Was the task easy or demanding, simple or complex?'];...
    ['작업 수행시 느낀 시간적 압박감의 정도',...
    newline,'How much time pressure did you feel due to the pace at which the tasks or task elements occurred? Was the pace slow or rapid?'];...
    ['작업 수행시 신체적 활동이 요구된 정도. 신체적으로 힘이 덜 들었나요, 많이 들었나요?',...
    newline, 'How much physical activity was required? Was the task easy or demanding, slack or strenuous?'];...
    ['작업을 수행하는 것이 성공적이었다고 느낀 정도',...
    newline, 'How successful were you in performing the task? How satisfied were you with your performance?'];...
    ['목표 달성을 위해 정신적, 신체적으로 노력한 정도',...
    newline,'How hard did you have to work (mentally and physically) to accomplish your level of performance?'];...
    ['작업 수행 중 짜증나고, 화가 나고, 스트레스를 받은 정도',...
    newline,'How irritated, stressed, and annoyed versus content, relaxed, and complacent did you feel during the task?']};



%%-- make recorder object
info = audiodevinfo;
recObj = audiorecorder(44100,16,1,1);
%% Pre (Control)
fprintf('Ready for Pre session (without video)\n')
fprintf('Run RDA_test_adaptive.exe and start stim without video\n')

pause;

MW = tcpclient('127.0.0.1',1800);
fopen(MW);

figure(Figs.h)
instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',54)

sess = 1;
set(Figs.h,'windowstate','minimized');
fwrite(MW,'51'); % trigger: session start
ClickCommand(1550,550); % move cursor
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

fprintf('session end\n')
pause;

ClickCommand(2500,100);
figure(Figs.h)
instruction(Figs.h,MW,Figs.ax,'이제 느낀 대로 응답해주세요. 스페이스바를  누르면 넘어갑니다.',52)

%-- NASA TLX
fwrite(MW,'31'); % trigger: NASA TLX response
response_rating = nasatlx_ratings(Figs.h2,loadtypes,descriptions);
TLX(sess,:) = response_rating;
set(Figs.h2,'WindowState','minimized')

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
set(Figs.h,'windowstate','minimized');
set(Figs.h2,'windowstate','minimized');
set(Figs.h3,'windowstate','minimized');

fprintf('Ready for Main session (with video)\n')
fprintf('Run RDA_test_adaptive.exe and start stim with video\n')



for sess = 2:Nsess-Npost

    fprintf(sprintf('Session%d (main %d)\n',sess,sess-1))
    instruction_hold(Figs.h3,[],Figs.ax_h3,'잠시 대기해주세요',50)

    MW = tcpclient('127.0.0.1',1800);
    fopen(MW);

    for iii = 1:3 %adjust cursor location
        ClickCommand(2500,100);
        MoveCommand(0,0);
    end
    ClickCommand(3800,500); % start stimuli
    fwrite(MW,'21'); % trigger: start stimuli

    figure(Figs.h)
    instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',54)
    set(Figs.h,'windowstate','minimized');

    record(recObj);
    fwrite(MW,'41'); % trigger: start record & session
    ClickCommand(1550,550); % move cursor

    for tr = 1:Ntr
        fprintf('\nSession%d Trial %d\n..',sess,tr)
        fprintf('Stim start in %ds\n',intervals(sess-1,tr))
        pause(intervals(sess-1,tr))
        target = targetlist(tr,sess);
        figure(Figs.h)
        targetpresent(Figs.h,MW,Figs.axAll,Figs.stim_img,target); % 3 sec

        %-- start BCI
        STARTBCI
        set(Figs.h,'windowstate','minimized');
        fwrite(MW,'57'); % trigger: start BCI
        fprintf('> BCI start\n')
        pause(18);
        fprintf('\nTrial %d end\n',tr)
    end

    fwrite(MW,'42'); % trigger: stop record & session
    stop(recObj);
    answer = getaudiodata(recObj);
    ansfile = [SubName,'_Sess',num2str(sess),'.wav'];
    audiowrite(ansfile, answer, recObj.SampleRate)

    fprintf('session end\n')
    pause;

    ClickCommand(2500,100);
    figure(Figs.h)
    instruction(Figs.h,MW,Figs.ax,'이제 느낀 대로 응답해주세요. 스페이스바를  누르면 넘어갑니다.',52)

    %-- NASA TLX
    fwrite(MW,'31'); % trigger: NASA TLX response
    response_rating = nasatlx_ratings(Figs.h2,loadtypes,descriptions);
    TLX(sess,:) = response_rating;
    set(Figs.h2,'WindowState','minimized')

    clear MW
end



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
set(Figs.h,'windowstate','minimized');
set(Figs.h2,'windowstate','minimized');
set(Figs.h3,'windowstate','minimized');
fprintf('Ready for Post session (without video)\n')
fprintf('Run RDA_test_adaptive.exe and start stim without video\n')

instruction_hold(Figs.h3,[],Figs.ax_h3,'잠시 대기해주세요',50)

MW = tcpclient('127.0.0.1',1800);
fopen(MW);

for iii = 1:3 %adjust cursor location
    ClickCommand(2500,100);
    MoveCommand(0,0);
end
ClickCommand(3800,500); % start stimuli
fwrite(MW,'21'); % trigger: start stimuli

for sess = Nsess-Npost+1:Nsess
    figure(Figs.h)
    instruction(Figs.h,MW,Figs.ax,'준비되면 스페이스바를 눌러주세요',54)
    set(Figs.h,'windowstate','minimized');
    fwrite(MW,'51'); % trigger: session start

    ClickCommand(1550,550); % move cursor
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

    fprintf('session end\n')
    pause;

    ClickCommand(2500,100);
    figure(Figs.h)
    instruction(Figs.h,MW,Figs.ax,'이제 느낀 대로 응답해주세요. 스페이스바를  누르면 넘어갑니다.',52)

    %-- NASA TLX
    fwrite(MW,'31'); % trigger: NASA TLX response
    response_rating = nasatlx_ratings(Figs.h2,loadtypes,descriptions);
    TLX(sess,:) = response_rating;
    set(Figs.h2,'WindowState','minimized')
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
%% TLX
global response_source

source_weight_comb = nasatlx_sourceofwl(Figs.h2,loadtypes,descriptions);
set(Figs.h2,'WindowState','minimized')
%-- calculate ratings
Nload = length(loadtypes);
sourceweight = zeros(Nload,1);
for n = 1:length(loadtypes)
    sourceweight(n) = sum(ismember(response_source,loadtypes{n}));
end
TLX(:,4) = 100 - TLX(:,4); % reverse PERFORMANCE rating
tlxresult.loadtype = loadtypes;
tlxresult.rating = TLX;
tlxresult.source = response_source;
tlxresult.source_weight_pair = source_weight_comb;
tlxresult.source_weight = sourceweight;
tlxresult.adjusted_rating = TLX.*sourceweight';
tlxresult.weighted_rating = sum(tlxresult.adjusted_rating,2)/15;

%% Finish
instruction_hold(Figs.h3,[],Figs.ax_h3,'실험이 모두 종료되었습니다',50)


vars.target = targetlist;
vars.intervals = intervals;
vars.subname = SubName;
vars.TLX = tlxresult;
save(SubName,'vars')

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
fprintf('==============\nAcc_pre  %.2f \n \nAcc_main1  %.2f \nAcc_main2  %.2f\nAcc_main3  %.2f\nAcc_main4  %.2f\nAcc_main5  %.2f\nAcc_main(all) %.2f\n\nAcc_post %.2f\n',...
    vars.Acc_pre,vars.Acc_main_sess(1),vars.Acc_main_sess(2),vars.Acc_main_sess(3),vars.Acc_main_sess(4),vars.Acc_main_sess(5),vars.Acc_mainAll,vars.Acc_post);

%-- save variables
save(SubName,'vars')
