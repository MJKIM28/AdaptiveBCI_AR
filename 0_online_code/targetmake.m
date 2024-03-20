%%
function axAll = targetmake(fig,winsize)
Side = winsize(4)/4;
axAll{1} = axes('Units','pixels','Position',[Side Side*2.5 Side Side],'color', [0 0 0],'xticklabel',{},'ticklength',[0 0],'parent',fig);
axAll{2} = axes('Units','pixels','Position',[winsize(3)-2*Side Side*2.5 Side Side],'color', [0 0 0],'xticklabel',{},'ticklength',[0 0],'parent',fig);
axAll{3} = axes('Units','pixels','Position',[Side Side*0.5 Side Side],'color', [0 0 0],'xticklabel',{},'ticklength',[0 0],'parent',fig);
axAll{4} = axes('Units','pixels','Position',[winsize(3)-2*Side Side*0.5 Side Side],'color', [0 0 0],'xticklabel',{},'ticklength',[0 0],'parent',fig);

set(axAll{1},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
set(axAll{2},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
set(axAll{3},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
set(axAll{4},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');

end
