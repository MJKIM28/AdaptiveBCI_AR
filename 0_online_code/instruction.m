function instruction(fig,com,ax,inst,trig1)
global clicked
clicked = 0;

set(fig,'windowstate','fullscreen','KeyPressFcn',@keypress);

th = text(ax,0.5,1,inst,'color',[1 1 1],'fontsize',45);
set(ax,'xticklabel',{},'ticklength',[0 0],'box','off','visible','off','xlim',[0 1],'ylim',[0 2]);

set(th,'visible','on','VerticalAlignment','middle','HorizontalAlignment', 'center')

if nargin > 5
    fwrite(com,num2str(trig1)); % trigger: instruction presented
end


while ~clicked
    pause(0.01);
end



% ylim([0 1])
% pause(duration);
cla(ax);

end

function keypress(src,event)
global clicked
if strcmp(event.Key,'space')
    clicked = 1;

end

end
