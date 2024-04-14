function topo_P300_partial(chdrawName,data,colors,conditions)
Ncon = length(conditions);

figure;
pb = [];
axes('Position',[0.0300    0.8    0.2134    0.1430])
hold on;
Fs = 500;
Times = -0.2:1/Fs:0.6-1/Fs;
drawaxis([-1.2 1.7],Times)
ylim([-1.2 1.7])
set(gca,'fontsize',12,'Box','off')
axis off

for ch = 1:length(chdrawName)
    if strcmp(chdrawName{ch},'Fz')
        axes('Position',[0.4108    0.7673    0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'FC1')
        axes('Position',[0.18 0.65  0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'FC2')
        axes('Position',[0.64 0.65  0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'Cz')
        axes('Position',[0.4108    0.55    0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'P3')
        axes('Position',[0.0300    0.3291    0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'Pz')
        axes('Position',[0.4108    0.3291    0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'P4')
        axes('Position',[ 0.78    0.3291    0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'O1')
        axes('Position',[ 0.18    0.1100    0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'Oz')
        axes('Position',[ 0.4108    0.1100    0.2134    0.1430])
    elseif strcmp(chdrawName{ch},'O2')
        axes('Position',[0.64    0.1100    0.2134    0.1430])
    end

    hold on

    if ismember('O',chdrawName{ch})
        ylimval = [-3.5 2];

    else
        ylimval = [-1.2 1.7];

    end
    ylim(ylimval)


    drawaxis(ylimval,Times)



    for c = 1:Ncon
        pb(c) = plot(Times,data(ch,:,c),'color',colors{c},'linewidth',2);


    end
    if ch == 1
        ylabel('Amplitude (\muV)')
        xlabel('Time (s)')
        legend(pb,conditions ...
            )
    else
        set(gca,'XTicklabel',[])
    end


    title(chdrawName{ch})

    set(gca,'fontsize',12,'Box','off')
    axis off
end
