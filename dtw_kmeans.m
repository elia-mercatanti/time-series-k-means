function [idx, new_C] = dtw_kmeans(X, k, init_C, error, max_iter)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% dtw_kmeans(X, k, initC, prec, maxIter)
%   Implementa l'algoritmo di clustering k-means utilizzando la distanza
%   DTW (Dynamic Time Wrap).
%
% Input:
% - X: Data Set.
% - k: numero di cluster da individuare.
% - init_C: matrice per l'inizializzazione dei centroidi. (Default random).
% - error: errore max con cui si vogliono ottenere i cluster. (Default 0)
% - max_iter: massimo numero di iterazioni a eseguire. (Default 100)
% Output:
% - new_idx: vettore di indici, indica il cluster di ogni oggetto.
% - new_C: matrice dei centroidi generati dall'algoritmo.

% verifica di alcune condizioni necessarie
[n, m] = size(X);
if nargin < 2
    error('Non sono stati inseriti tutti gli input.');
elseif (~mod(k,1) == 0) || (k == 1) || (k > n-1)
    error("Il numero di cluster deve rispettare l'intervallo [2, n-1]");
elseif exist('init_C', 'var') && (size(init_C,1)~=k || size(init_C,2)~=m)
    error('La matrice dei centroidi iniziali deve avere dimensioni:k x m');
end

% assegnazione dei parametri di default
if ~exist('max_iter', 'var')
    max_iter = 100;
end
if ~exist('error', 'var')
    error = 0;
end
if ~exist('init_C', 'var')
    old_C = zeros(k, m);
    ts_index = randperm(n, k);
    for i = 1:k
        old_C(i,:) = X(ts_index(i), :);
    end   
else
    old_C = init_C;
end

% inizializzo le variabili iniziali
distance = zeros(1, k);
ts_to_C = zeros(1, n);
new_C = zeros(k, m);
idx = zeros(n, 1);
convergence = false;
iter = 0;

% loop principale
while ~convergence
    % inizializzo alcune matrici per il ricalcolo dei cluster
    clusters_sum = zeros(k, m);
    clusters_num = zeros(k, 1);
    
    % assegno ogni serie temporale al centroide piu' vicino
    for i = 1:n
        for j = 1:k
            % distanza fra serie temporale e centroide
            distance(j) = dtw(X(i,:), old_C(j,:));
        end
        % ottengo l'indice della distanza minima (centroide piu' vicino)
        [min_dist, min_idx] = min(distance);
        % assegno la serie temporale al centroide piu' vicino
        idx(i) = min_idx;
        % ottengo la distanza della s. t. dal suo centroide
        ts_to_C(i) = min_dist;
        
        % mi calcolo via via la somma delle s. t. di uno stesso cluster
        clusters_sum(min_idx, :) = clusters_sum(min_idx, :) + X(i,:);
        clusters_num(min_idx) = clusters_num(min_idx) + 1;
    end
    
    % ricalcolo i centroidi facendo al media delle s. t. assegnate ad essi
    for i = 1:k
        % controllo se un cluster e' rimasto vuoto.
        if clusters_num(i) == 0
            % ricavo l'indice della s. t. piu' distante dal suo centroide 
            % e la assegno al cluster vuoto
            [~, farthest_ts] = max(ts_to_C);
            original_c = idx(farthest_ts);
            idx(farthest_ts) = i;
            % aggiornando le altre variabili di conseguenza
            new_C(i, :) = X(farthest_ts, :);
            clusters_num(i) = 1;
            clusters_num(original_c) = clusters_num(original_c) - 1;
            clusters_sum(i, :) = X(farthest_ts, :);
            clusters_sum(original_c, :) = clusters_sum(original_c, :) ...
                                          - X(farthest_ts, :);
            new_C(original_c, :) = clusters_sum(original_c, :) / ...
                                   clusters_num(original_c);
               
            % mi ricalcolo il vettore delle distanze delle dai centroidi
            ts_to_C(farthest_ts) = dtw(X(farthest_ts,:), new_C(i,:));
            original_clust_ts = find(idx == original_c);
            for j = 1:length(original_clust_ts)
                ts_idx = original_clust_ts(j);
                ts_to_C(ts_idx) = dtw(X(ts_idx, :), new_C(original_c, :));
            end
        else
            new_C(i, :) = clusters_sum(i, :) / clusters_num(i);
        end
    end
    iter = iter + 1;
    
    % controllo criterio convergenza
    if iter == max_iter | (all(abs(old_C - new_C)) <= error) %#ok<OR2>
        convergence = true;
    else
        old_C = new_C;
    end
end
end