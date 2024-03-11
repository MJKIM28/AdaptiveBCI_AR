
function topo_P300(ch_locs,tar_data,ntar_data,c,ylimit,LegendName)
% target data shape = (channel, time points)
% nontarget data shape = (nstim, channel, time points)
channel = ch_locs;
x = [-100:299] / 500;
% ylimit = [-2 2];

% figure()
if any(strcmp(channel,'FP1')) == 1
    subplot(7,7,3)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FP1')),:)), 'Color',c  , 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FP1')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('Fp1')
end
if any(strcmp(channel,'FPZ')) == 1
    subplot(7,7,4)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FPZ')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FPZ')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('Fpz')
end
if any(strcmp(channel,'FP2')) == 1
    subplot(7,7,5)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FP2')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FP2')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('Fp2')
end
if any(strcmp(channel,'F7')) == 1
    subplot(7,7,9)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'F7')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'F7')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('F7')
end
if any(strcmp(channel,'F3')) == 1
    subplot(7,7,10)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'F3')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'F3')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('F3')
end
if any(strcmp(channel,'FZ')) == 1
    subplot(7,7,11)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FZ')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FZ')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('Fz')
end
if any(strcmp(channel,'F4')) == 1
    subplot(7,7,12)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'F4')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'F4')),:))', 'k-')
    hold on
    plot([0 0],ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('F4')
end
if any(strcmp(channel,'F8')) == 1
    subplot(7,7,13)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'F8')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'F8')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('F8')
end
if any(strcmp(channel,'FT9')) == 1
    subplot(7,7,15)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FT9')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FT9')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('FT9')
end
if any(strcmp(channel,'FC5')) == 1
    subplot(7,7,16)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FC5')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FC5')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('FC5')
end
if any(strcmp(channel,'FC1')) == 1
    subplot(7,7,17)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FC1')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FC1')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('FC1')
end
if any(strcmp(channel,'FC2')) == 1
    subplot(7,7,19)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FC2')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FC2')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('FC2')
end
if any(strcmp(channel,'FC6')) == 1
    subplot(7,7,20)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FC6')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FC6')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit)
        xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('FC6')
end
if any(strcmp(channel,'FT10')) == 1
    subplot(7,7,21)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'FT10')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'FT10')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('FT10')
end
if any(strcmp(channel,'T7')) == 1
    subplot(7,7,23)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'T7')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'T7')),:))', 'k-')
    hold on
    plot([0 0],ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('T7')
end
if any(strcmp(channel,'C3')) == 1
    subplot(7,7,24)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'C3')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'C3')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('C3') 
end
if any(strcmp(channel,'CZ')) == 1
    subplot(7,7,25)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'CZ')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'CZ')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('Cz')
end
if any(strcmp(channel,'C4')) == 1
    subplot(7,7,26)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'C4')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'C4')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('C4')
end
if any(strcmp(channel,'T8')) == 1
    subplot(7,7,27)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'T8')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'T8')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('T8')
end
if any(strcmp(channel,'CP5')) == 1
    subplot(7,7,30)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'CP5')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'CP5')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('CP5')
end
if any(strcmp(channel,'CP1')) == 1
    subplot(7,7,31)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'CP1')),:)),'Color',c,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'CP1')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('CP1')
end
if any(strcmp(channel,'CP2')) == 1
    subplot(7,7,33)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'CP2')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'CP2')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('CP2')
end
if any(strcmp(channel,'CP6')) == 1
    subplot(7,7,34)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'CP6')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'CP6')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('CP6')
end
if any(strcmp(channel,'P7')) == 1
    subplot(7,7,37)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'P7')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'P7')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('P7')
end
if any(strcmp(channel,'P3')) == 1
    subplot(7,7,38)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'P3')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'P3')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('P3')
end
if any(strcmp(channel,'PZ')) == 1
    subplot(7,7,39)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'PZ')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'PZ')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('Pz')
end
if any(strcmp(channel,'P4')) == 1
    subplot(7,7,40)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'P4')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'P4')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('P4')
end
if any(strcmp(channel,'P8')) == 1
    subplot(7,7,41)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'P8')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'P8')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('P8')
end
if any(strcmp(channel,'O1')) == 1
    subplot(7,7,45)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'O1')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'O1')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('O1')
end
if any(strcmp(channel,'OZ')) == 1
    subplot(7,7,46)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'OZ')),:)),'Color',c, 'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'OZ')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('Oz')
end
if any(strcmp(channel,'O2')) == 1
    subplot(7,7,47)
    plot(x, squeeze(tar_data(find(strcmp(channel, 'O2')),:)),'Color',c ,'LineWidth', 1.5)
    hold on
    plot(x, squeeze(ntar_data(:,find(strcmp(channel, 'O2')),:))', 'k-')
    hold on
    plot([0 0], ylimit,'k','LineWidth', 2) % onset
    ylim(ylimit);    xlim([min(x),max(x)]);
    xlabel('sec(s)')
    title('O2')
end

lgd=legend(LegendName,'Location','Southwest');
lgd.FontSize = 20;
lgd.Color = 'none';

fig = subplot(7,7,48);
plot([0 0], ylimit,'Color','#999999','LineStyle',':','LineWidth', 2) % onset
box off
xlabel('sec(s)')
ylabel('Amplitude(mV)')
fig.LineWidth = 1.5;
 ylim(ylimit)
 xlim([-0.2 0.6])
set(fig, 'color', 'none')
end