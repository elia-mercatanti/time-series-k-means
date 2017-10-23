% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% Test: 
%   Script per testare l'efficenza delll'algoritmo 'TSkmeans' rispetto ad 
%   altri algoritmi di clustering per serie temporali.

% algoritmi disponibili
algorithm = {'Time Series K-Means', 'Euclid K-Means', 'Pearson K-Means',...
   'City Block K-Means', 'Cosine K-Means', 'DTW K-Means', 'DTW K-Medoids', 'Euclid K-Medoids'}; 
             
% ottengo quali algoritmi voglio testare tra quelli disponibili
[chosen_alg, algorithm_num] = get_algorithm(algorithm);

% ricavo il Data Set su cui testare gli algoritmi e il vettore delle classi
[X, true_idx, set_name] = get_data_set();
[n, m] = size(X);

% ricavo il numero di cluster da ricercare
k = input('Inserire il numero di cluster da ricercare: ');

% calcolo il global scatter del Data Set e il parametro alpha ('TSkmeans')
if chosen_alg(1) == 1
    parameter = input('Inserire il valore del parametro (alpha/gs): ');
    [gs] = global_scatter(X);
    alpha = gs * parameter;
end

% ottengo il num. di iterazioni del test
iter = input('Inserire il numero di iterazioni del test: ');

% inizializzo la matrice dei tempi di esecuzione e le matrici degli indici 
% interni/esterni per la valutazione dei cluster generati
[timings, sse_indices, ssb_indices, sil_indices, pur_indices, ...
  fs_indices, ran_indices, nmi_indices] = deal(zeros(algorithm_num, iter));

% inizializzo gli array di celle dei centroidi e dei pesi generati ad ogni
% iterazione, la matrice dei centroidi iniziali e gli indici dei cluster
init_C = zeros(k, m);
idx = zeros(n, algorithm_num);
C = cell(1, algorithm_num);
W = cell(1, algorithm_num);
W(2:algorithm_num) = {0};

% loop principale per il test dei vari algoritmi di clustering
fprintf('\n##################### Inizio dei Test #####################\n');
fprintf('Iterazioni eseguite:   0');
for i = 1:iter
    % genero casualmente i centroidi iniziali per l'i-esima iter. del test
    ts_index = randperm(n, k);
    for j = 1:k
        init_C(j,:) = X(ts_index(j),:);
    end
    
    % eseguo gli algoritmi selezionati per il test
    for j = 1:algorithm_num
        alg_tested = chosen_alg(j);
        tic
        switch alg_tested
            case 1
                [U, clust, W{j}] = ts_kmeans(X, k, alpha, init_C);
                
                % calcolo il vettore idx corretto per 'TSkmeans'
                [row, col] = find(U == 1);
                clust_idx = sortrows([row, col], 1);
                clust_idx = clust_idx(:, 2);
            case 2
                [clust_idx, clust] = kmeans(X, k, 'Start', init_C);
            case 3
                [clust_idx, clust] = kmeans(X, k, 'Distance', ...
                                           'correlation', 'Start', init_C);
            case 4
                [clust_idx, clust] = kmeans(X, k, 'Distance', ...
                                             'cityblock', 'Start', init_C);
            case 5
                [clust_idx, clust] = kmeans(X, k, 'Distance', 'cosine', ...
                                                          'Start', init_C);                     
            case 6
                [clust_idx, clust] = dtw_kmeans(X, k, init_C);
            case 7
                [clust_idx, clust] = kmedoids(X, k, 'Distance', ...
                                               @dtw_dist, 'Start', init_C);
            case 8
                [clust_idx, clust] = kmedoids(X, k, 'Start', init_C);                               
        end
        timings(j, i) = toc;
        idx(:,j) = clust_idx;
        C{j} = clust;
        
        % calcolo varie metriche/indici per la valutazione degli algoritmi
        sse_indices(j, i) = sse(X, clust, clust_idx);
        ssb_indices(j, i) = ssb(X, clust, clust_idx);
        sil_indices(j, i) = mean(silhouette(X, clust_idx));
        pur_indices(j, i) = purity(true_idx, clust_idx);
        fs_indices(j, i) = f_score(true_idx, clust_idx);
        ran_indices(j, i) = rand_index(true_idx, clust_idx);
        nmi_indices(j, i) = nmi(true_idx, clust_idx);
    end
    fprintf('\b\b\b%3d', i);
end
fprintf('\n');
fprintf('##################### Test Completati #####################\n\n');

fprintf('################### Risultati dei Test ####################\n\n');
% calcolo il valore medio delle metriche/indici e dei tempi di esecuzione
timings = mean(timings, 2);
sse_indices = mean(sse_indices, 2);
ssb_indices = mean(ssb_indices, 2);
tss_metrics = sse_indices + ssb_indices;
sil_indices = mean(sil_indices, 2);
pur_indices = mean(pur_indices, 2);
fs_indices = mean(fs_indices, 2);
ran_indices = mean(ran_indices, 2);
nmi_indices = mean(nmi_indices, 2);

% creo la tabella dei tempi di esec., la salvo in un file e la mostro
fprintf(['## Tempi di Esecuzione (media sulle %u iterazioni) degli ' ...
                    'algoritmi sul Data Set "%s" ##\n\n'], iter, set_name);
timings_table = table(algorithm(chosen_alg)', timings, 'VariableNames', ...
                                                 {'Algorithm', 'Timing'});
disp(timings_table);
fprintf('## Indici dei risultati per il Data Set "%s" ##\n', set_name);

% creo la tabella degli indici interni, la salvo in un file e la mostro
fprintf('\nIndici Interni (media sulle %u iterazioni):\n\n', iter);
var_names = {'Algorithm', 'SSE', 'SSB', 'TSS', 'Silhouette'};
inter_idx_table = table(algorithm(chosen_alg)', sse_indices, ...
        ssb_indices, tss_metrics, sil_indices, 'VariableNames', var_names);
disp(inter_idx_table);

% creo la tabella degli indici esterni, la salvo in un file e la mostro
fprintf('Indici Esterni (media sulle %u iterazioni):\n\n', iter);
var_names = {'Algorithm', 'Purity', 'F_Score', 'Rand_Index', 'NMI'};
exter_idx_table = table(algorithm(chosen_alg)', pur_indices, fs_indices,... 
                    ran_indices, nmi_indices, 'VariableNames', var_names);
disp(exter_idx_table);

% decido se voglio salvare i risultati in file excel
answer = input('Vuoi salvare i risultati?(Y o N): ', 's');
if answer == 'Y'
    writetable(timings_table, ['Risultati\timings-' set_name ' (alpha='...
                                                 num2str(alpha) ').xlsx']);
    writetable(inter_idx_table, ['Risultati\internal_idx-' set_name ...
                                       '(alpha=' num2str(alpha) ').xlsx']);
    writetable(exter_idx_table, ['Risultati\external_idx-' set_name ...
                                       '(alpha=' num2str(alpha) ').xlsx']);
end

% decido se voglio disegnare i grafici dei risultati dell'ultima iterazione
answer = input(['Mostrare i grafici dei risultati dell''ultima ' ...
                                            'iterazione?(Y o N): '], 's');
if answer == 'Y'
    % dimensioni matrice dei grafici per la visualizzazione dei cluster
    request = ('Dimensioni del grafico dei clusters ([righe, colonne]): ');
    clust_plot = input(request);
    if length(clust_plot) ~= 2
        error('Deve essere rispettata la sintassi: [righe, colonne].');
    end
    % disegno i grafici dei vari algoritmi
    for i = 1:algorithm_num
        plotTitle = ['Risultati del "' algorithm{chosen_alg(i)} '" sul '...
                                               'Data Set "' set_name '"'];
        plot_results(X, C{i}, W{i}, idx(:, i), clust_plot, chosen_alg(i)...
                                                             , plotTitle);
    end
end

function [chosen_alg, algorithm_num] = get_algorithm(algorithm)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% get_algorithm(algorithm): 
%   Funzione per chiedere all'utente quali saranno gli algoritmi, fra 
%   quelli disponibili, sui quali effettuare i test.
%
% Input:
% - algorithm: array di celle, contiene i nomi degli algoritmi disponibili.
% Output:
% - chosen_alg: array binario, contiene gli algoritmi scelti dall'utente.
% - algorithm_num: numero totale degli algoritmi scelti.

algorithm_num = length(algorithm);
fprintf('Algoritmi disponibili: ');
for i = 1:algorithm_num
    fprintf('%u) %s, ', i, algorithm{i});
end
fprintf('\b\b.\n\n');
chosen_alg = input('Indicare su quali algoritmi eseguire il test: ');
if size(chosen_alg) ~= algorithm_num
    error('Rispettare la sintassi: [alg1{1,0},alg2{1,0},...,algNum{1,0}]');
end
chosen_alg = find(chosen_alg);
algorithm_num = length(chosen_alg);
end

function [X, true_idx, set_name] = get_data_set()
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% get_data_set(): 
%   Funzione per ricavare il Data Set scelto dall'utente e per ottenere il
%   vettore di indici che indica la vera classe di ogni serie temporale.
%
% Input:
% - algorithm: array di celle, contiene i nomi degli algoritmi disponibili.
% Output:
% - X: matrice delle serie temporali.
% - true_idx: vettore di indici, indica la vera classe di ogni s. t..
% - set_name: array binario, contiene gli algoritmi scelti dall'utente.

% chiedo all'utente di inserire il nome del Data Set che si vuole testare
request = ('Inserire il nome del Data Set su cui testare gli algoritmi:');
set_name = input([request ' '], 's');
dir_path = ['UCR_Data_Sets\' set_name];
data_path = [dir_path '\' set_name '_'];

% cerso il Data Set e unisco il suo Test Set con il relativo Training Set
if exist(dir_path, 'dir') == 7
    X = load([data_path  'TEST']);
    X = [X; load([data_path 'TRAIN'])];                                          
elseif exist(set_name, 'file') == 2
    X = load(set_name);
else
    error('Il Data Set inserito non esiste.');
end

% ricavo il vettore che mi indica la classe di appartenenza di ciascuna 
% serie temporale (UCR) e lo separo dal data set vero e proprio.
true_idx = X(:, 1);
if any(ismember(0, true_idx))
    true_idx =  true_idx + 1;
end
X = X(:, 2:end);
end

function plot_results(X, C, W, idx, clust_plot, algorithm_id, plot_title)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% plot_results(X, C, W, idx, clust_plot, algorithm_id, plot_title): 
%   Disegna i grafici dei risultati degli algoritmi di clustering.
%
% Input:
% - X: matrice delle serie temporali.
% - C: matrice dei centroidi di ogni cluster.
% - W: matrice dei pesi di ogni time stamp dei cluster.
% - idx: vettore di indici, indica il cluster di ogni serie temporale.
% - clust_plot: dimensione del grafico dei clusters.
% - algorithm_id: identifica a quale algoritmo corrispondono i grafici.
% - plot_title: contiene il titolo dei grafici.

% ricavo alcuni parametri fissi utilizzati per il disegno dei grafici
[n, m] = size(X);
k = size(C, 1);
legend_info = cellstr(num2str((1:k)', 'Cluster n.%d'));

% disegno il grafico delle serie temporali
limits = plot_data_set(X, n, m, plot_title);

% disegno il grafico dei centroidi generati
y_lim = plot_centroids(C, legend_info, algorithm_id, limits);

% disegno il grafico dei cluster trovati, con le serie che vi appartengono
plot_clusters(X, C, n, m, k, idx, algorithm_id, limits, y_lim)

% disegno il grafico dei pesi di ogni time stamps per ogni cluster
if W ~= 0
    plot_weights(W, m, legend_info);
end

% disegno il grafico di ogni cluster trovato, in modo separato
plot_isolate_clusters(X, C, m, k, idx, clust_plot, plot_title);
end

function [limits] = plot_data_set(X, n, m, plot_title)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% plot_data_set(X, n, m, plot_title): 
%   Disegna il grafico delle serie temporali passate in input.
%
% Input:
% - X: matrice delle serie temporali.
% - n: numero di serie temporali.
% - m: numero totale di time stamps per ogni serie temporale.
% - plot_title: contiene il titolo dei grafici.
% Output:
% - limits: limiti rispetto a x ed y del grafico.

figure('Name', plot_title);
set(gcf, 'units', 'normalized', 'outerposition', [0 0 1 1]);
subplot(2,2,1);
plot(X');
grid on
xlim([1 m]);
title(['Serie Temporali Esaminate, N=', num2str(n)]);
xlabel('Time Stamp');
ylabel('Valore');
limits = axis;
end

function [yl] = plot_centroids(C, legend_info, algorithm_id, limits)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% plot_centroids(C, legend_info, algorithm_id, limits): 
%   Disegna il grafico dei centroidi generati.
%
% Input:
% - C: matrice dei centroidi di ogni cluster.
% - legend_info: stringhe per generare la legenda del grafico.
% - algorithm_id: identifica a quale algoritmo corrispondono i grafici.
% - limits: limiti rispetto a x ed y del grafico.
% Output:
% - yl: limite rispetto alle sole y del grafico.

yl = zeros(1,2);
subplot(2,2,2);
plot(C', 'LineWidth', 1.5);
grid on
title('Centroidi Generati');
xlabel('Time Stamp');
ylabel('Valore');
legend(legend_info);
if algorithm_id ~= 3
    axis(limits);
else
    y_lim = ylim;
    if y_lim(1) <= limits(3)
        yl(1) = y_lim(1);
    else
        yl(1) = limits(3);
    end
    if y_lim(2) >= limits(4)
        yl(2) = y_lim(2);
    else
        yl(2) = limits(4);
    end
    xlim(limits(1:2));
    ylim([yl(1) yl(2)]);
end
end

function plot_clusters(X, C, n, m, k, idx, algorithm_id, limits, y_lim)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% plot_clusters(X, C, n, m, k, idx, algorithm_id, limits): 
%   Disegna il grafico dei cluster con le serie che vi appartengono.
%
% Input:
% - X: matrice delle serie temporali.
% - C: matrice dei centroidi di ogni cluster.
% - n: numero di serie temporali.
% - m: numero totale di time stamps per ogni serie temporale.
% - k: numero totale di cluster.
% - idx: vettore di indici, indica il cluster di ogni serie temporale.
% - algorithm_id: identifica a quale algoritmo corrispondono i grafici.
% - limits: limiti rispetto a x ed y del grafico.
% - y_lim: limite rispetto alle sole y del grafico.

subplot(2,2,3);
hold on
color_vector = hsv(k);
for i = 1:n
    plot(1:m, X(i,:), 'col', color_vector(idx(i),:));
end
for i = 1:k
    plot(1:m, C(i,:), 'k-x', 'LineWidth', 1.5);
end
grid on
title('Grafico di Appartenenza ai Cluster');
xlabel('Time Stamp');
ylabel('Valore');
hold off
if algorithm_id ~= 3
    axis(limits);
else
    xlim(limits(1:2));
    ylim(y_lim);
end
end

function plot_weights(W, m, legend_info)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% plot_weights(W, m, legend_info): 
%   Disegna il grafico dei pesi di ogni time stamps per ogni cluster.
%
% Input:
% - W: matrice dei pesi di ogni time stamp dei cluster.
% - m: numero totale di time stamps per ogni serie temporale.
% - legend_info: stringhe per generare la legenda del grafico.

subplot(2,2,4);
plot(W', 'LineWidth', 1.5, 'DisplayName', 'Cluster');
grid on
xlim([1 m]);
title('Pesi Generati per ogni Cluster');
xlabel('Time Stamp');
ylabel('Peso');
legend(legend_info);
end

function plot_isolate_clusters(X, C, m, k, idx, clust_plot, plot_title)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% plot_isolate_clusters(X, C, m, k, clust_plot, plot_title): 
%   Disegna il grafico di ogni cluster trovato, in modo separato.
%
% Input:
% - X: matrice delle serie temporali.
% - C: matrice dei centroidi di ogni cluster.
% - m: numero totale di time stamps per ogni serie temporale.
% - k: numero totale di cluster.
% - idx: vettore di indici, indica il cluster di ogni serie temporale.
% - clust_plot: dimensione del grafico dei clusters.
% - plot_title: contiene il titolo dei grafici.

figure('Name', plot_title);
set(gcf,'units', 'normalized', 'outerposition',[0 0 1 1]);
for i = 1:k
    subplot(clust_plot(1), clust_plot(2), i);
    hold on
    ts_idx = find(idx == i);
    for j = 1:length(ts_idx)
        plot(1:m, X(ts_idx(j), :));
    end
    plot(1:m, C(i,:), 'k-x', 'LineWidth', 1.5);
    grid on
    xlim([1 m]);
    title(['Cluster n.', num2str(i)]);
    xlabel('Time Stamp');
    ylabel('Valore');
    hold off
end
end

function d2 = dtw_dist(XI, XJ)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% dtw_dist(XI, XJ): 
%   Funzione utiizzata per specificare la distanza 'DTW' nel k-medoids.

m = size(XJ, 1);
result = zeros(m, 1);
for k = 1:m
    result(k) = dtw(XI, XJ(k,:));
end
d2 = result;
end