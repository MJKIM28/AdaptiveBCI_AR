% h = figure;
% screensize =  get(0, 'ScreenSize');
% set(h,'color','k','menubar','none','toolbar','none','position',screensize)

function tlxresult = NASATLX(h)
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;

%-- give cursor to participant
mouse.mouseMove(2500, 500)


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


set(h,'windowstate','fullscreen');

%% plot questions _ ratings
% global response_rating

response_rating = nasatlx_ratings(h,loadtypes,descriptions);

%% plot questions _ source of workload
global response_source

source_weight_comb = nasatlx_sourceofwl(h,loadtypes,descriptions);

%-- give cursor back to experimenter
mouse.mouseMove(500, 500)

%% calculation - weight, adjusted rating, weighted rating
Nload = length(loadtypes);
sourceweight = zeros(Nload,1);
for n = 1:length(loadtypes)
   sourceweight(n) = sum(ismember(response_source,loadtypes{n}));
end
response_rating(4,:) = 100 - response_rating(4,:); % reverse PERFORMANCE rating
tlxresult.loadtype = loadtypes;
tlxresult.rating = response_rating;
tlxresult.source = response_source;
tlxresult.source_weight_pair = source_weight_comb;
tlxresult.source_weight = sourceweight;
tlxresult.adjusted_rating = response_rating.*sourceweight;
tlxresult.weighted_rating = sum(tlxresult.adjusted_rating)/15;
end