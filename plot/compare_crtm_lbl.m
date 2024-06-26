clf;
addpath('./');

subplot(1,2,1);
%crtm
data = load('crtm_wf.n18_hirs4.with_o3');
plot(data(:,2), data(:,1),'marker','o','markeredgecolor','b','color','b');
set(gca,'yscale','log','ydir','reverse','ylim',[0.1 1000], 'xlim', [-0.1 0.5]);

%lbl
data2 = load('/Users/cda/Desktop/WeiFuncProj/lbl_gasabs/src/fort.1000');
hold on;
plot(data2(:,2),data2(:,1),'marker','o','markeredgecolor','r','color','r');
set(gca,'yscale','log','ydir','reverse','ylim',[0.1 1000], 'xlim', [-0.1 0.5]);

subplot(1,2,2);
%crtm
data = load('crtm_wf.n18_hirs4.with_o3');
plot(data(:,2)./max(data(:,2)), data(:,1),'marker','o','markeredgecolor','b','color','b');
set(gca,'yscale','log','ydir','reverse','ylim',[0.1 1000], 'xlim', [-0.1 1.1]);

%lbl
data2 = load('/Users/cda/Desktop/WeiFuncProj/lbl_gasabs/src/fort.1000');
hold on;
plot(data2(:,2)./max(data2(:,2)),data2(:,1),'marker','o','markeredgecolor','r','color','r');
set(gca,'yscale','log','ydir','reverse','ylim',[0.1 1000], 'xlim', [-0.1 1.1]);
