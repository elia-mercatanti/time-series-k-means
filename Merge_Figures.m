% Load saved figures
f1=hgload('Grafici\FaceFour-Alpha.fig');
f2=hgload('Grafici\ECGFiveDays-Alpha.fig');
f3=hgload('Grafici\ItalyPowerDemand-Alpha.fig');
f4=hgload('Grafici\synthetic_control-Alpha.fig');
f5=hgload('Grafici\ArrowHead-Alpha.fig');
alpha_gs = [0.05, 0.2:0.2:1, 2:2:10, 50:50:500, 1000];
val_num = length(alpha_gs);

% Prepare subplots
figure
h(1)=subplot(2,3,1);
grid on
title('(A) FaceFour');
ylabel('Misure');
xlim([log(alpha_gs(1)) log(alpha_gs(val_num))]);

h(2)=subplot(2,3,2);
grid on
title('(B) ECGFiveDays');
xlim([log(alpha_gs(1)) log(alpha_gs(val_num))]);

h(3)=subplot(2,3,3);
grid on
title('(C) ItalyPowerDemand');
xlabel('Ln(\alpha/gs)');
xlim([log(alpha_gs(1)) log(alpha_gs(val_num))]);

h(4)=subplot(2,3,4);
grid on
title('(D) SynteticControl');
xlabel('Ln(\alpha/gs)');
ylabel('Misure');
xlim([log(alpha_gs(1)) log(alpha_gs(val_num))]);

h(5)=subplot(2,3,5);
grid on
title('(E) ArrowHead');
xlabel('Ln(\alpha/gs)');
xlim([log(alpha_gs(1)) log(alpha_gs(val_num))]);

% Paste figures on the subplots
copyobj(allchild(get(f1,'CurrentAxes')),h(1));
copyobj(allchild(get(f2,'CurrentAxes')),h(2));
copyobj(allchild(get(f3,'CurrentAxes')),h(3));
copyobj(allchild(get(f4,'CurrentAxes')),h(4));
copyobj(allchild(get(f5,'CurrentAxes')),h(5));

% Add legends
legend('Purity', 'F-Score', 'Rand Index', 'NMI', 'Location', 'best');