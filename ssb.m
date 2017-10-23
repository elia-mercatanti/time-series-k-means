function [result] = ssb(X, C, idx)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
% 
% ssb(X, C, idx):
%      Calcola il Sum of Squares Between clusters (SSB).
%
% Input:
% - X: matrice delle serie temporali.
% - C: array di celle per i centroidi dei cluster.
% - idx: vettore che indica l'appartenenza di ogni s. t. al suo cluster.
% Output:
% - result: Sum of Squares Between clusters (SSB).

result = 0;
k = size(C,1);
overallC = sum(X)/size(X,1);
for i = 1:k
    dimC = length(find(idx == i));
    result = result + dimC*norm(C(i,:) - overallC)^2;
end
end