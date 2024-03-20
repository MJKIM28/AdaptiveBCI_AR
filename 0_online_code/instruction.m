function instruction(fig,com,ax,inst,duration,trig)

set(fig,'windowstate','fullscreen');

th = text(ax,0.5,1,inst,'color',[1 1 1],'fontsize',45);
set(ax,'xticklabel',{},'ticklength',[0 0],'box','off','visible','off','xlim',[0 1],'ylim',[0 2]);

set(th,'visible','on','VerticalAlignment','middle','HorizontalAlignment', 'center')

if nargin > 5
    fwrite(com,num2str(trig)); % trigger: instruction
end

% ylim([0 1])
pause(duration);
cla(ax);

end
