% h = figure;
% screensize =  get(0, 'ScreenSize');
% set(h,'color','k','menubar','none','toolbar','none','position',screensize)

function response_rating = nasatlx_ratings(h,types,descriptions)
% global response_rating
% global axesloc
% global objmade

Ntype = length(types);
axeslist = cell(Ntype,1);
x1 = 0:10:100;
Nx = length(x1);
y1_1 = zeros(1,Nx);
y1_2 = ones(1,Nx);
y1_2(round(Nx/2)) = 1.5;
x2 = 5:10:100;
Nx2 = length(x2);
y2_1 = zeros(1,Nx2);
y2_2 = 0.5*ones(1,Nx2);

response_rating = zeros(Ntype,1);
axesloc = zeros(Ntype,1);
objmade = cell(Ntype,1);
figure(h);
for n = 1:Ntype
    
    axeslist{n} = subplot(Ntype,1,n,'Color','k','ButtonDownFcn',@nasatlx_ButtonDown_ClickRate);
    axesloc(n) = axeslist{n}.Position(2); % y value of position
    
    plot(axeslist{n},[x1; x1],[y1_1; y1_2],'w','ButtonDownFcn',@nasatlx_ButtonDown_ClickRate);
    hold on;
    plot(axeslist{n},x1,y1_1,'w','ButtonDownFcn',@nasatlx_ButtonDown_ClickRate);
    plot(axeslist{n},[x2; x2],[y2_1; y2_2],'w','ButtonDownFcn',@nasatlx_ButtonDown_ClickRate);
    axis(axeslist{n},'off');
    if ismember('Performance',types{n})
         text(0,2,[descriptions{n}] ,'color','w','fontsize',15,'HitTest','off');       
        text(0,0,'Poor','HorizontalAlignment','center','VerticalAlignment','top','color','w','HitTest','off');
        text(100,0,'Good','HorizontalAlignment','center','VerticalAlignment','top','color','w','HitTest','off');
    else
        text(0,2,descriptions{n},'color','w','fontsize',15,'HitTest','off');
        text(0,0,'Low','HorizontalAlignment','center','VerticalAlignment','top','color','w','HitTest','off');
        text(100,0,'High','HorizontalAlignment','center','VerticalAlignment','top','color','w','HitTest','off');
    end
    ylim([0 3])
    title(axeslist{n},types{n},'color','w','fontsize',20,'HitTest','off');
    
end

nasatlx_nextpage(h,axeslist);


%%
function nasatlx_ButtonDown_ClickRate(hObject, eventdata)
% global axesloc
% global response_rating
% global objmade

    aa = get(gca,'CurrentPoint');
    score = 5*round(aa(1)/5);
    aaa = get(gca);
    axesselected = find(axesloc == aaa.Position(2)); % find which question answered
    
    if ~isempty(response_rating(axesselected)) % if this question already answered 
        delete(objmade{axesselected}); % erase previous answer      
    end

    response_rating(axesselected) = score;  
    objmade{axesselected} = plot(gca,score,0.5,'rv','Markersize',5,'linewidth',5); % update the answer
    

end

end
