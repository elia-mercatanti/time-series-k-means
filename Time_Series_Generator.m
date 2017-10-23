time_stamps = 20;
time_series = 5;
legend_info = cellstr(num2str((1:time_series)', 'Serie Temporale n.%d'));

x = 1:time_stamps;
y = zeros(time_series, time_stamps);
for i = 1:time_series
    y(i,:) = cumsum(randn(length(x),1));
end

plot(x, y, 'LineWidth', 2);
xlim([1 time_stamps]);
legend(legend_info)
xlabel('Time Stamp');
ylabel('Valore');
grid on;