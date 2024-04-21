function instruction_hold(fig,com,ax,inst,trig1)

set(fig,'windowstate','fullscreen');

th = text(ax,0.5,1,inst,'color',[1 1 1],'fontsize',45);
set(ax,'xticklabel',{},'ticklength',[0 0],'box','off','visible','off','xlim',[0 1],'ylim',[0 2]);

set(th,'visible','on','VerticalAlignment','middle','HorizontalAlignment', 'center')

if nargin > 5
    fwrite(com,num2str(trig1)); % trigger: instruction presented
end

fprintf('press anything if ready:\n')
pause;

cla(ax);
set(fig,'WindowState','minimized')


end
