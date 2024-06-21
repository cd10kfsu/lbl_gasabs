cd /Users/cda/Desktop/WeiFuncProj/src;
lev = load('fort.100');
lay = load('fort.200');

clf;
%temperature
subplot(1,3,1);
hold on;
plot(lev(:,2),lev(:,1)/100,'marker','o');
plot(lay(:,2),lay(:,1)/100,'marker','o','markeredgecolor','r','color','r');

%H2O
subplot(1,3,2);
hold on;
plot(lev(:,4),lev(:,1)/100,'marker','o');
plot(lay(:,4),lay(:,1)/100,'marker','o','markeredgecolor','r','color','r');


%CO2
subplot(1,3,3);
hold on;
plot(lev(:,5),lev(:,1)/100,'marker','o');
plot(lay(:,5),lay(:,1)/100,'marker','o','markeredgecolor','r','color','r');
