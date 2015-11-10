addpath('./');
data = load('fort.1000');

clf;
subplot(1,2,1);
hold on;
plot(data(:,2)*1/max(data(:,2)),data(:,1),'marker','o','markeredgecolor','r','color','r');
set(gca,'yscale','log','ydir','reverse','ylim',[1 1000], 'xlim', [-0.1 1.1]);

