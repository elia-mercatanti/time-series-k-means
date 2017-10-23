function [result] = sse(X, C, idx)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
% 
% see(X, C, idx):
%      Calcola il Sum of Squarred Error (SSE).
%
% Input:
% - X: matrice delle serie temporali.
% - C: array di celle per i centroidi dei cluster.
% - idx: vettore che indica l'appartenenza di ogni s. t. al suo cluster.
% Output:
% - result: Sum of Squarred Error (SSE).

result = 0;
k = size(C,1);
for i = 1:k
    tsIdx = find(idx == i);
    for j = 1:length(tsIdx)
        result = result + norm(C(i,:) - X(tsIdx(j),:))^2;
    end
end
end