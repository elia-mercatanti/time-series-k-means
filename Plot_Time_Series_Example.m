data_set = load('time_series_example');
data_set = data_set(:, 2:end);
legend_info = cellstr(num2str((1:size(data_set,1))', 'Serie Temporale n.%d'));

figure('Name', 'Data set in esame');
plot(data_set', 'LineWidth', 2);
grid on;
legend(legend_info);
title('Serie Temporali del set in esame');
xlabel('Time Stamp');
ylabel('Valore');
xlim([1 size(data_set,2)]);