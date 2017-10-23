function [rand] = rand_index(true_idx, idx)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% rand_index(true_idx, idx): 
%   Calcola la metrica 'Rand-Index', che valuta l'accuratezza dei cluster.
% Reference: http://refhub.elsevier.com/S0020-0255(16)30379-6/sbref0021
%
% Input:
% - true_idx: vettore di indici, indica la vera classe di ogni oggetto.
% - idx: vettore di indici, indica il cluster di ogni oggetto.
% Output:
% - rand: valore finale di della metrica 'Rand-Index'.

% inizializzo alcuni parametri iniziali
n = length(idx);
a = 0;
b = 0;

% calcolo il Rand-Index
for i = 1:n-1
    for j = i+1:n
        if (idx(i) == idx(j)) && (true_idx(i) == true_idx(j))
            a = a + 1;
        end
        if (idx(i) ~= idx(j)) && (true_idx(i) ~= true_idx(j))
            b = b + 1;
        end
    end
end
rand = (a + b) / (n*(n-1)/2);
end