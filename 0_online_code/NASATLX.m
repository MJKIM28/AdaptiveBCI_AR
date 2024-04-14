% h = figure;
% screensize =  get(0, 'ScreenSize');
% set(h,'color','k','menubar','none','toolbar','none','position',screensize)

function tlxresult = NASATLX(h)
import java.awt.Robot;
import java.awt.event.*;
mouse = Robot;

%-- give cursor to participant
mouse.mouseMove(2500, 500)


loadtypes = {'������ ���̵� (Mental Demand)';'�ð��� ���̵� (Temporal Demand)';'��ü�� ���̵� (Physical Demand)';...
    '���밨 (Performance)';'��� ���� (Effort)';'������ (Frustration)'};
descriptions = {['�۾� ����� ������/������ Ȱ���� �䱸�� ����. �۾��� ��������, ���������? �ܼ��߳���, �����߳���?',...
    newline,'How much mental and perceptual activity was required? Was the task easy or demanding, simple or complex?'];...
        ['�۾� ����� ���� �ð��� �йڰ��� ����',...
        newline,'How much time pressure did you feel due to the pace at which the tasks or task elements occurred? Was the pace slow or rapid?'];...
    ['�۾� ����� ��ü�� Ȱ���� �䱸�� ����. ��ü������ ���� �� �������, ���� �������?',...
    newline, 'How much physical activity was required? Was the task easy or demanding, slack or strenuous?'];...
    ['�۾��� �����ϴ� ���� �������̾��ٰ� ���� ����',...
    newline, 'How successful were you in performing the task? How satisfied were you with your performance?'];...
    ['��ǥ �޼��� ���� ������, ��ü������ ����� ����',...
    newline,'How hard did you have to work (mentally and physically) to accomplish your level of performance?'];...
    ['�۾� ���� �� ¥������, ȭ�� ����, ��Ʈ������ ���� ����',...
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