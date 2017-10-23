load('synthetic_data_set');
data_set = sortrows(synthetic_data_set);
data_set = data_set(:, 2:end);
Clust1 = data_set(1:100, 1:end);
Clust2 = data_set(101:200, 1:end);
Clust3 = data_set(201:300, 1:end);

figure('Name', 'Data set Sintetico');
subplot(2,2,1);
plot(data_set');
grid on
title('(A) Serie Temporali del Set Sintetico');
xlabel('Time Stamp');
ylabel('Valore');
xlim([1 15]);

subplot(2,2,2);
plot(Clust1');
grid on
title('(B) Serie Temporali del Cluster 1');
xlabel('Time Stamp');
ylabel('Valore');
xlim([1 15]);

subplot(2,2,3);
plot(Clust2');
grid on
title('(C) Serie Temporali del Cluster 2');
xlabel('Time Stamp');
ylabel('Valore');
xlim([1 15]);

subplot(2,2,4);
plot(Clust3');
grid on
title('(D) Serie Temporali del Cluster 3');
xlabel('Time Stamp');
ylabel('Valore');
xlim([1 15]);