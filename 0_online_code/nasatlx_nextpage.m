function nasatlx_nextpage(h,ax)
global clicked
clicked = 0;
A = axes(h,'position',[0.9 0.45 0.01 0.01],'color','k');
hold on
text(A,0,0,'NEXT >>','color','w','fontweight','bold','fontsize',20,'verticalalignment','bottom','ButtonDownFcn',@nasatlx_ButtonDown_Pressnext);
axis(A,'off');

while ~clicked
    pause(0.01);
end


    function nasatlx_ButtonDown_Pressnext(hObject, eventdata)
        % clear the page: delete axes object
        for n = 1:length(ax)
        delete(ax{n}); 
        end
        delete(A); 
        clicked = 1;
    end
end

