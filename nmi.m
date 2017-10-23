function [nmi] = nmi(true_idx, idx)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% nmi(true_idx, idx): 
%   Calcola la Normalized Mutual Information (NMI), una metrica per la 
%   valutazione della qualita' dei cluster.
% Reference: http://refhub.elsevier.com/S0020-0255(16)30379-6/sbref0021
%
% Input:
% - true_idx: vettore di indici, indica la vera classe di ogni oggetto.
% - idx: vettore di indici, indica il cluster di ogni oggetto.
% Output:
% - nmi: valore finale della metrica 'NMI'.

% inizializzo alcuni parametri iniziali e ottengo la matrice di confusione
n = length(idx);
mutual_info = 0;
clusters_entropy = 0;
classes_entropy = 0;
confusion = confusionmat(idx, true_idx);
k = size(confusion, 1);

% calcolo la NMI
np = sum(confusion, 2);
nq = sum(confusion);
for p = 1:k
    for q = 1:k
        npq = confusion(p, q);
        if npq ~= 0
            mutual_info = mutual_info + npq*log2((n*npq)/(np(p)*nq(q)));
        end        
    end
    clusters_entropy = clusters_entropy + np(p)*log2(np(p)/n);
    classes_entropy = classes_entropy + nq(p)*log2(nq(p)/n);
end
nmi = mutual_info/sqrt(clusters_entropy*classes_entropy);
end