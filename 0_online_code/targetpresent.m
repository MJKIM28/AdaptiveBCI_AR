%%
function targetpresent(fig,com,ax,stim,target)

set(fig,'windowstate','fullscreen');
fprintf('> target....%d\n',target);



% fill(ax{1},[0 0 1 1],[0 1 1 0],[98 88 254]/255);
% fill(ax{2},[0 0 1 1],[0 1 1 0],[98 88 254]/255);
% fill(ax{3},[0 0 1 1],[0 1 1 0],[98 88 254]/255);
% fill(ax{4},[0 0 1 1],[0 1 1 0],[98 88 254]/255);

set(ax{1},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
set(ax{2},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
set(ax{3},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
set(ax{4},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');

fwrite(com,num2str(target));
switch target
  
 case 1
        imshow(stim{1,2},'InitialMagnification','fit','Parent',ax{1})
        set(ax{1},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        
        imshow(stim{2,1},'InitialMagnification','fit','Parent',ax{2})
        set(ax{2},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{3,1},'InitialMagnification','fit','Parent',ax{3})
        set(ax{3},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{4,1},'InitialMagnification','fit','Parent',ax{4})
        set(ax{4},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        
    case 2
        imshow(stim{2,2},'InitialMagnification','fit','Parent',ax{2});
        set(ax{2},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        
        imshow(stim{1,1},'InitialMagnification','fit','Parent',ax{1})
        set(ax{1},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{3,1},'InitialMagnification','fit','Parent',ax{3})
        set(ax{3},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{4,1},'InitialMagnification','fit','Parent',ax{4})
        set(ax{4},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');        
    case 3
        imshow(stim{3,2},'InitialMagnification','fit','Parent',ax{3});
        set(ax{3},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        
        imshow(stim{1,1},'InitialMagnification','fit','Parent',ax{1})
        set(ax{1},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{2,1},'InitialMagnification','fit','Parent',ax{2})
        set(ax{2},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{4,1},'InitialMagnification','fit','Parent',ax{4})
        set(ax{4},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        
    case 4
        imshow(stim{4,2},'InitialMagnification','fit','Parent',ax{4});
        set(ax{4},'xticklabel',{},'ticklength',[0 0],'yticklabel',{},'box','off','visible','off');

        imshow(stim{1,1},'InitialMagnification','fit','Parent',ax{1})
        set(ax{1},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{2,1},'InitialMagnification','fit','Parent',ax{2})
        set(ax{2},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
        imshow(stim{3,1},'InitialMagnification','fit','Parent',ax{3})
        set(ax{3},'xticklabel',{},'ticklength',[0 0],'box','off','visible','off');
end


pause(3);

cla(ax{1});
cla(ax{2});
cla(ax{3});
cla(ax{4});

set(fig,'windowstate','minimized');
% fwrite(com,'56'); % trigger: target presentation
end

