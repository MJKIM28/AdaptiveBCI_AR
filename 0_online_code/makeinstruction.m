%%
function ax = makeinstruction(fig,winsize)
Side = winsize(4)/4;
ax = axes('Units','pixels','Position',[Side Side winsize(3)-2*Side winsize(4)-2*Side],'color', [0 0 0],'xticklabel',{},'ticklength',[0 0],'parent',fig);

set(ax,'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
end
