function [fs] = f_score(true_idx, idx)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% f_score(true_idx, idx): 
%   Calcola la metrica 'F-Score', ovvero, una combinazione delle metriche 
%   'Precision' e 'Recall'.
% Reference: http://refhub.elsevier.com/S0950-7051(14)00266-4/h0015
%
% Input:
% - true_idx: vettore di indici, indica la vera classe di ogni oggetto.
% - idx: vettore di indici, indica il cluster di ogni oggetto.
% Output:
% - fs: valore finale della metrica 'F-Score'.

% inizializzo alcuni parametri iniziali e ottengo la matrice di confusione
n = length(idx);
fs = 0;
confusion = confusionmat(true_idx, idx);
k = size(confusion, 1);
temp_fs = zeros(1, k);

% calcolo la F-Score
nj = sum(confusion, 2);
nl = sum(confusion);
for j = 1:k
    for l = 1:k
        njl = confusion(j, l);
        temp_fs(l) = (2*(njl/nj(j))*(njl/nl(l)))/((njl/nj(j))+(njl/nl(l)));
    end
    fs = fs + (nj(j)/n)*max(temp_fs);
end
end