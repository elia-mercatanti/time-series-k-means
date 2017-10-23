function [pur] = purity(true_idx, idx)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% purity(true_idx, idx): 
%   Calcola la metrica 'Purity', ovvero, la percentuale degli oggetti 
%   correttamente clusterizzati.
% Reference: https://goo.gl/SWN3De
%
% Input:
% - true_idx: vettore di indici, indica la vera classe di ogni oggetto.
% - idx: vettore di indici, indica il cluster di ogni oggetto.
% Output:
% - pur: valore finale della metrica 'Purity'.

% inizializzo alcuni parametri iniziali e ottengo la matrice di confusione
n = length(idx);
confusion = confusionmat(true_idx, idx);

% calcolo la Purity
pur = sum(max(confusion)) / n;
end


