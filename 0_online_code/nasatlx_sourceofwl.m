% h = figure;
% screensize =  get(0, 'ScreenSize');
% set(h,'color','k','menubar','none','toolbar','none','position',screensize)


function combination = nasatlx_sourceofwl(h,loadtypes,descriptions)
global response_source
global answer
global clicked

Ntype = length(loadtypes);

combination = nchoosek(1:Ntype,2);
Ncomb = size(combination,1);

id = randperm(Ncomb);
combination = combination(id,:);

response_source = cell(Ncomb,1);
figure(h)
clicked = 0;
for n = 1:Ncomb
pair = combination(n,:);

subplot(2,2,[1,2],'color','k');
text(0.5,0.7,[loadtypes{pair(1)} ,newline, descriptions{pair(1)}],'color','w','horizontalalignment','center','fontsize',15);
text(0.5,0.4,[loadtypes{pair(2)} ,newline,descriptions{pair(2)}],'color','w','horizontalalignment','center','fontsize',15);
text(0.5,0,'아래의 두 요인 중 작업 수행시 느낀 작업 부하에 더 크게 영향을 끼친 요인을 클릭하세요.','color','w','horizontalalignment','center','fontsize',25);
axis off 

subplot(2,2,3,'color','k');
text(0.5,0.5,loadtypes{pair(1)},'color','w','ButtonDownFcn',@nasatlx_ButtonDown_choose,'horizontalalignment','center','fontsize',30);
axis off

subplot(2,2,4,'color','k');
text(0.5,0.5,loadtypes{pair(2)},'color','w','ButtonDownFcn',@nasatlx_ButtonDown_choose,'horizontalalignment','center','fontsize',30);
axis off

while ~clicked
    pause(0.001); % pause until the response is given
end
clicked = 0;
response_source{n} = answer;
clf;
end
% end