

load(['E:\1_연구\분석코드\2024_adaptive_BCI\Data\Pilot\Train\Epoch\Subtest20.mat'])

colors = [0 128 128; 255 105 180;0 191 255]./255;
Fs = 500;
Times = -0.2:1/Fs:0.6-1/Fs;
ind = Times>0;
figure; plot(Times(ind),Epoch.dat(ind,22,1,1)+15,'color',colors(1,:),'LineWidth',2)
hold on; plot(Times(ind),Epoch.dat(ind,23,1,1),'color',colors(2,:),'LineWidth',2)
hold on; plot(Times(ind),Epoch.dat(ind,24,1,1)-15,'color',colors(3,:),'LineWidth',2)
set(gca,'ytick',[],'Box','off','FontSize',13)
xlabel("Time (ms)")