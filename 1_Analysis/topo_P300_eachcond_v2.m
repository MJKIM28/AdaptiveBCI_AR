
function topo_P300_eachcond_v2(ch_locs,data,colors,conditions)
% data shape = (channel, time points, Ncondition)

Ncon = size(data,3);
Ch =                                 {'Fp1'    'Fpz'    'Fp2'    'F7'    'F3'    'Fz'    'F4'    'F8' ...
    'FT9'    'FC5'    'FC1'    'FC2'    'FC6'    'FT10'    'T7'    'C3'...
    'Cz'    'C4'    'T8'    'CP5'    'CP1'    'CP2'    'CP6'    'P7'...
    'P3'    'Pz'    'P4'    'P8'    'O1'    'Oz'    'O2'};
ChP300 = { 'C3' 'Cz'    'C4'       'CP1'    'CP2'   ...
    'P3'    'Pz'    'P4'     'O1'    'Oz'    'O2'};
channel = ch_locs;
x = [-99:300] / 500;
ylimit = [-2 2];

fig=figure();
for ich = 1: size(Ch,2)
    if any(strcmp(channel(ich),'Fp1')) == 1
        subplot(7,7,3)
    elseif any(strcmp(channel(ich),'Fpz')) == 1
        subplot(7,7,4)
    elseif any(strcmp(channel(ich),'Fp2')) == 1
        subplot(7,7,5)
    elseif any(strcmp(channel(ich),'F7')) == 1
        subplot(7,7,9)
    elseif any(strcmp(channel(ich),'F3')) == 1
        subplot(7,7,10)
    elseif any(strcmp(channel(ich),'Fz')) == 1
        subplot(7,7,11)
    elseif any(strcmp(channel(ich),'F4')) == 1
        subplot(7,7,12)
    elseif any(strcmp(channel(ich),'F8')) == 1
        subplot(7,7,13)
    elseif any(strcmp(channel(ich),'FT9')) == 1
        subplot(7,7,15)
    elseif any(strcmp(channel(ich),'FC5')) == 1
        subplot(7,7,16)
    elseif any(strcmp(channel(ich),'FC1')) == 1
        subplot(7,7,17)
    elseif any(strcmp(channel(ich),'FC2')) == 1
        subplot(7,7,19)
    elseif any(strcmp(channel(ich),'FC6')) == 1
        subplot(7,7,20)
    elseif any(strcmp(channel(ich),'FT10')) == 1
        subplot(7,7,21)
    elseif any(strcmp(channel(ich),'T7')) == 1
        subplot(7,7,23)
    elseif any(strcmp(channel(ich),'C3')) == 1
        subplot(7,7,24)
    elseif any(strcmp(channel(ich),'Cz')) == 1
        subplot(7,7,25)
    elseif any(strcmp(channel(ich),'C4')) == 1
        subplot(7,7,26)
    elseif any(strcmp(channel(ich),'T8')) == 1
        subplot(7,7,27)
    elseif any(strcmp(channel(ich),'CP5')) == 1
        subplot(7,7,30)
    elseif any(strcmp(channel(ich),'CP1')) == 1
        subplot(7,7,31)
    elseif any(strcmp(channel(ich),'CP2')) == 1
        subplot(7,7,33)
    elseif any(strcmp(channel(ich),'CP6')) == 1
        subplot(7,7,34)
    elseif any(strcmp(channel(ich),'P7')) == 1
        subplot(7,7,37)
    elseif any(strcmp(channel(ich),'P3')) == 1
        subplot(7,7,38)
    elseif any(strcmp(channel(ich),'Pz')) == 1
        subplot(7,7,39)
    elseif any(strcmp(channel(ich),'P4')) == 1
        subplot(7,7,40)
    elseif any(strcmp(channel(ich),'P8')) == 1
        subplot(7,7,41)
    elseif any(strcmp(channel(ich),'O1')) == 1
        subplot(7,7,45)
    elseif any(strcmp(channel(ich),'Oz')) == 1
        subplot(7,7,46)
    elseif any(strcmp(channel(ich),'O2')) == 1
        subplot(7,7,47)
    end
    hold on;
%     if ismember(Ch(ich),ChP300)
%         A =fill([0.25 0.25 0.5 0.5],[ylimit fliplr(ylimit)],'k','facealpha',.2,'linestyle','none');
%     end
b  = [];
for con = 1:Ncon
    b(con) = plot(x, squeeze(data(find(strcmp(channel, Ch(ich))),:,con)),'Color',colors{con}, 'LineWidth', 2); hold on;
end

    
    str = '#999999';
    color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
    C = plot([0 0], ylimit,'Color',color','LineStyle',':','LineWidth', 1); % onset
    ylim(ylimit)
    set(gca,'XTickLabel',[],'YTickLabel',[])
    title(strrep(Ch(ich),'Z','z'), 'FontSize', 16)
    box on
end

lgd=legend(b,conditions,'Location','Southwest');
lgd.FontSize = 20;
lgd.Color = 'none';

fig = subplot(7,7,48);
hold on;
plot([0 0], ylimit,'Color',color','LineStyle',':','LineWidth', 2);

% fill([0.25 0.25 0.5 0.5],[ylimit fliplr(ylimit)],'k','facealpha',.2,'linestyle','none');

box off
xlabel('Time(s)')
ylabel('Amplitude(\muV)')
fig.LineWidth = 1.5;
ylim(ylimit)
xlim([-0.2 0.6])

set(fig, 'color', 'none')

end