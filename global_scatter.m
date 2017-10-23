function [gs] = global_scatter(X)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% global_scatter(X):
%   Funzione che ha il compito di calcolare il global scatter della
%   matrice delle serie temporali.
% Reference: http://dx.doi.org/10.1016/j.ins.2016.05.040
%
% Input:
% - X: matrice delle serie temporali.
% Output:
% - gs: stima finale del global scatter.

[n, m] = size(X);
gs = 0;
for i = 1:n
    for j = 1:m
        z = sum(X(:,j))/n;
        gs = gs + (X(i,j) - z)^2;
    end
end
end

