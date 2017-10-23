% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% Test_Alpha: 
%   Script per testare le prestazioni dell'algoritmo TSkmeans al variare
%   del parametro alpha. Consente di scegliere il corretto valore di alpha.

% algoritmi disponibili
algorithm = {'Time Series K-Means'}; 

% ricavo il Data Set su cui testare gli algoritmi e il vettore delle classi
[X, true_idx, set_name] = get_data_set();
[n, m] = size(X);

% ricavo il numero di cluster da ricercare
k = input('Inserire il numero di cluster da ricercare: ');

% calcolo il global scatter del Data Set e il parametro alpha ('TSkmeans')
[gs] = global_scatter(X);

% ottengo il num. di iterazioni del test
iter = input('Inserire il numero di iterazioni del test: ');

% inizializzo le matrici degli indici esterni per la valutazione dei 
% cluster generati
[pur_indices, fs_indices, ran_indices, nmi_indices] = deal(zeros(1, iter));

% inizializzo i valori di (alpha/gs) per i test da effetuare, da -2.9957 
% (alpha/gs=10^3) a 6.9078 (alpha/gs=10^3), e i valori del parametro alpha
% fare riferimento all'articolo per chiarimenti
alpha_gs = [0.05, 0.2:0.2:1, 2:2:10, 50:50:500, 1000];
alpha = alpha_gs*gs;

% inizializzo le variabili utilizzate per salvare le metriche risultanti
% dalle varie iterazioni.
val_num = length(alpha_gs);
measures = zeros(4, val_num);

% loop principale per il test del parametro alpha
fprintf('\n##################### Inizio dei Test #####################\n');
fprintf('Valori testati:   0');
for i = 1:val_num  
    for j = 1:iter  
        % eseguo l'algoritmo 'TSkmeans' con il particolare valore di alpha
        [U, ~, ~] = ts_kmeans(X, k, alpha(i));

        % calcolo il vettore idx corretto per 'TSkmeans'
        [row, col] = find(U == 1);
        clust_idx = sortrows([row, col], 1);
        clust_idx = clust_idx(:, 2);

        % calcolo varie metriche/indici per la valutazione degli algoritmi
        pur_indices(j) = purity(true_idx, clust_idx);
        fs_indices(j) = f_score(true_idx, clust_idx);
        ran_indices(j) = rand_index(true_idx, clust_idx);
        nmi_indices(j) = nmi(true_idx, clust_idx);   
    end
    
    % calcolo il valore medio delle metriche/indici
    measures(1, i) = mean(pur_indices, 2);
    measures(2, i) = mean(fs_indices, 2);
    measures(3, i) = mean(ran_indices, 2);
    measures(4, i) = mean(nmi_indices, 2);
    fprintf('\b\b\b%3d', i);
end
fprintf('\n');
fprintf('##################### Test Completati #####################\n\n');

% disegno il grafico degli indici ottenuti al variare del valore di alpha
figure('Name', ['Prestazioni del TSkmeans al variare del parametro ' ...
                                '\alpha per il Data Set "' set_name '"']);
plot(log(alpha_gs), measures, 'LineWidth', 2, 'Marker', 'x');
grid on
xlabel('Ln(\alpha/gs)');
ylabel('Misure');
xlim([log(alpha_gs(1)) log(alpha_gs(val_num))]);
legend('Purity', 'F-Score', 'Rand Index', 'NMI', 'Location', 'best');

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
request = (['Inserire il nome del Data Set su cui trovare il valore ' ...
                                        'ottimale del parametro alpha:']);
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