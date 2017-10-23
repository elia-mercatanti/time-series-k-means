% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
%
% Synthetic_Data_Generator: 
%   Script per generare un data set sintetico con cluster ben definiti
%   (simile all'esempio dell'articolo).

% inizializzo i cluster e il data set
C1 = zeros(100, 15);
C2 = zeros(100, 15);
C3 = zeros(100, 15);
X = zeros(300, 15);
membership = zeros(300, 1);

% genero le serie temporali per ogni cluster
for i = 1:100
    C1(i, 1) = (0.55-0.4).*rand + 0.4;
    C1(i, 2) = (0.8-0.6).*rand + 0.6;
    C1(i, 3) = (0.9-0.7).*rand + 0.7;
    C1(i, 4) = (0.6-0.45).*rand + 0.45;
    C1(i, 5) = (0.75-0.6).*rand + 0.6;
    C1(i, 6:15) = rand(1, 10);
    C2(i, 1:5) = rand(1, 5);
    C2(i, 6) = (0.3-0.1).*rand + 0.1;
    C2(i, 7) = (0.4-0.2).*rand + 0.2;
    C2(i, 8) = (0.45-0.3).*rand + 0.3;
    C2(i, 9) = (0.3-0.2).*rand + 0.2;
    C2(i, 10) = (0.4-0.3).*rand + 0.3;
    C2(i, 11:15) = rand(1, 5);
    C3(i, 1:10) = rand(1, 10);
    C3(i, 11) = (1-0.85).*rand + 0.85;
    C3(i, 12) = (0.85-0.65).*rand + 0.65;
    C3(i, 13) = (1-0.80).*rand + 0.80;
    C3(i, 14) = (0.85-0.7).*rand + 0.7;
    C3(i, 15) = (0.95-0.825).*rand + 0.825;
end

% creo il data set ordinato per cluster
X(1:100, 1:15) = C1;
X(101:200, 1:15) = C2;
X(201:300, 1:15) = C3;

% vettore di etichette di classe
membership(1:100) = 1;
membership(101:200) = 2;
membership(201:300) = 3;

% riordino casualmente le righe e ottengo il set sintetico finale X
X = [membership X];
X = X(randperm(end), :);
save('synthetic_data_set', 'X', '-ascii');
