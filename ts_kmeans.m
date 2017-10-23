function [U, Z, W] = ts_kmeans(X, k, alpha, init_Z)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
% 
% ts_kmeans(X, k, alpha, init_Z):
%   Implementa l'algoritmo di smooth subspace clustering "TSkmeans" per 
%   serie temporali.
% Reference: http://dx.doi.org/10.1016/j.ins.2016.05.040
%
% Input:
% - X: matrice di n serie temporali X={X_1,X_2,...X_n), rappresentate da
%      vettori. Ogni vettore di una serie temporale X_i={x_i1,x_i2,...x_im)
%      e' caratterizzato da m valori che coincidono con m time stamps.
% - k: numero di cluster da individuare.
% - alpha: parametro alpha usato per bilanciare gli effeti tra la
%          dispersione degli oggetti all'interno dei cluster e tra la
%          levigatezza dei pesi dei time stamps.
% - init_Z: vettore utilizzato per inizializzare i centroidi.
% Output:
% - U: matrice di membership binaria n x k, dove l'elemento u_ip=1 indica
%      che una serie temporale i e' stata assegnata al cluster p,
%      altrimenti non e' assegnata a p.
% - Z: matrice dei centroidi dei k cluster, rappresentate da
%      vettori Z={Z_1,Z_2,...Z_k) che contengono le serie temporali che
%      meglio rappresentano l’andamento delle serie appartenenti ad un
%      determinato cluster, i centroidi sono quindi serie temporali.
% - W: W={W_1,W_2,...W_k) è un set di k vettori che rappresentatno i pesi
%      dei vari time stamps in ogni cluster. Il valore dell'elemento w_pj
%      indica il peso del j-esimo time stamp per il p-esimo cluster.

[n, m] = size(X);

% verifica di alcune condizioni necessarie
if nargin < 3
    error('Non sono stati inseriti tutti gli input necessari.');
elseif (~mod(k,1) == 0) || (k == 1) || (k > n-1)
    error("Il numero di cluster deve rispettare l'intervallo [2, n-1]");
elseif exist('init_Z', 'var') && (size(init_Z,1)~=k || size(init_Z,2)~= m)
    error('La matrice dei centroidi iniziali deve avere dimensioni:k x m');
end

% inizializzazione
[Z, W, H, Aeq, beq, lb, ub] = initialize(X, k, alpha, n, m);
if exist('init_Z', 'var')
    Z = init_Z;
end
U0 = zeros(n, k);
convergence = false;

% loop principale
while ~convergence
    % ricavo la matrice U
    [U] = solve_U(X, Z, W, k, n);
    
    % ricavo la matrice Z
    [Z] = solve_Z(X, U, k, m);
    
    % ricavo la matrice W
    [W] = solve_W(X, U, Z, k, n, m, H, Aeq, beq, lb, ub);
    
    % controllo il criterio di convergenza
    if isequal(U0, U)
        convergence = true;
    else
        U0 = U;
    end
end
end

function [Z, W, H, Aeq, beq, lb, ub] = initialize(X, k, alpha, n, m)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
% 
% initialize(X, k, alpha, n, m):
%   Inizializza le matrici iniziali Z0 e W0, che verranno utilizzate nell'
%   algoritmo di clustering. Inoltre, ha il compito di calcolare la matrice
%   hessiana H e i vincoli del problema di programmazione quadratica.
% Reference: http://dx.doi.org/10.1016/j.ins.2016.05.040
%
% Input:
% - X: matrice di n serie temporali X={X_1,X_2,...X_n), rappresentate da
%      vettori. Ogni vettore di una serie temporale X_i={x_i1,x_i2,...x_im)
%      e' caratterizzato da m valori che coincidono con m time stamps.
% - k: numero di cluster da individuare.
% - alpha: parametro alpha usato per bilanciare gli effeti tra la
%          dispersione degli oggetti all'interno dei cluster e tra la
%          levigatezza dei pesi dei time stamps.
% - n: numero di serie temporali.
% - m: numero totale di time stamps per ogni serie temporale.
% Output:
% - Z: matrice dei centroidi dei k cluster, rappresentate da
%      vettori Z0={Z_1,Z_2,...Z_k) che contengono le serie temporali che
%      meglio rappresentano l’andamento delle serie appartenenti ad un
%      determinato cluster, i centroidi sono quindi serie temporali.
% - W: W0={W_1,W_2,...W_k) è un set di k vettori che rappresentatno i pesi
%      dei vari time stamps in ogni cluster. Il valore dell'elemento w_pj
%      indica il peso del j-esimo time stamp per il p-esimo cluster.
% - H:  matrice simmetrica che rappresenta la forma quadratica di
%       1/2*w'*H*w+f'*w., utilizzata dal tool 'quadprog'.
% - Aeq: matrice che rappresenta i coefficenti lineari del vincolo
%        Aeq*w=beq.
% - beq: rappresenta il vettore costante nel vincolo Aeq*x = beq.
% - lb: vettore che rappresenta il limite inferiore nel vincolo lb<=w<=ub.
% - ub: vettore che rappresenta il limite superiore nel vincolo lb<=w<=ub.

% inizializzo le matrici principali
Z = zeros(k, m);
W = zeros(k, m);

% vengono generati k indici di serie temporali in X in modo random
seriesIndex = randperm(n, k);
for i = 1:k
    % le serie puntate dagli indici vengono scelte come centroidi iniziali.
    Z(i,:) = X(seriesIndex(i),:);
    
    % genero random i pesi iniziali di ogni time stamps per ogni cluster
    W(i,:) = randfixedsum(m, 1, 1, 0, 1)';
end

% ricavo la matrice simmetrica hessiana H per il tool quadprog
H = zeros(k*m, k*m);
for p = 1:k
    for j = 1:m-1
        index1 = j+m*(p-1);
        index2 = j+1+m*(p-1);
        H(index1,index1) = H(index1,index1) + 1;
        H(index2,index2) = H(index2,index2) + 1;
        H(index1,index2) = H(index1,index2) - 1;
        H(index2,index1) = H(index2,index1) - 1;
    end
end
H = H * alpha;

% ricavo il primo vincolo da applicare al problema di programmazione
% quadratica (Sommatoria da j=1 a m di w_pj deve essere 1)
Aeq = zeros(k, k*m); 
for i = 1:k
    Aeq(i, (i-1)*m+1:i*m) = 1;
end
beq = ones(k, 1);

% ricavo il secondo vincolo da applicare al problema di programmazione
% quadratica (lb = 0 <= w_pj <= 1 = ub)
lb = zeros(k*m, 1);
ub = ones(k*m, 1);
end

function [U] = solve_U(X, Z, W, k, n)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
% 
% solve_U(X, Z, W, k, n, m):
%   Ricava la matrice U minimizzando la funzione obiettivo dell'algoritmo 
%   di clustering prendendo come costanti le matrici W e Z e utilizzando 
%   la formula (3).
% Reference: http://dx.doi.org/10.1016/j.ins.2016.05.040
%
% Input:
% - X: matrice di n serie temporali X={X_1,X_2,...X_n), rappresentate da
%      vettori. Ogni vettore di una serie temporale X_i={x_i1,x_i2,...x_im)
%      e' caratterizzato da m valori che coincidono con m time stamps.
% - Z: matrice dei centroidi dei k cluster, rappresentate da
%      vettori Z={Z_1,Z_2,...Z_k) che contengono le serie temporali che
%      meglio rappresentano l’andamento delle serie appartenenti ad un
%      determinato cluster, i centroidi sono quindi serie temporali.
% - W: W={W_1,W_2,...W_k) è un set di k vettori che rappresentatno i pesi
%      dei vari time stamps in ogni cluster. Il valore dell'elemento w_pj
%      indica il peso del j-esimo time stamp per il p-esimo cluster.
% - k: numero di cluster da individuare.
% - n: numero di serie temporali.
% Output:
% - U: matrice di appartenenza binaria n x k, dove l'elemento u_ip=1 indica
%      che una serie temporale i e' stata assegnata al cluster p,
%      altrimenti non e' assegnata a p.

% fissato Z e W trovo la matrice U utilizzando la formula (3)
U = zeros(n, k);
distance = zeros(1, k);
dist_to_Z = zeros(1, n);
for i = 1:n
    for p = 1:k
        distance(p) = sum(W(p, :).*((X(i, :)-Z(p, :)).^2));
    end
    [min_dist, min_dist_idx] = min(distance);
    U(i, min_dist_idx) = 1;
    dist_to_Z(i) = min_dist;
end
        
% controllo se un cluster e' rimasto vuoto, nel caso ci inserisco una s. t.
empty_clusters = find(~any(U)); 
if ~isempty(empty_clusters)
    for i = 1:length(empty_clusters)
        % viene scelta la s. t. piu' lontana dal suo centroide
        [~, farthest_ts] = max(dist_to_Z);
        original_cluster = find(U(farthest_ts, :));
        
        % aggiorno le principali variabili di conseguenza
        U(farthest_ts, empty_clusters(i)) = 1;
        U(farthest_ts, original_cluster) = 0; %#ok<FNDSB>
        Z(empty_clusters(i), :) = X(farthest_ts, :);
        dist_to_Z(farthest_ts) = 0;
    end
end
end

function [Z] = solve_Z(X, U, k, m)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
% 
% solve_Z(X, U, k, m):
%   Ricava la matrice Z minimizzando la funzione obiettivo dell'algoritmo 
%   di clustering prendendo come costanti le matrici U e W e utilizzabdo
%   la formula (5).
% Reference: http://dx.doi.org/10.1016/j.ins.2016.05.040
%
% Input:
% - X: matrice di n serie temporali X={X_1,X_2,...X_n), rappresentate da
%      vettori. Ogni vettore di una serie temporale X_i={x_i1,x_i2,...x_im)
%      e' caratterizzato da m valori che coincidono con m time stamps.
% - U: matrice di appartenenza binaria n x k, dove l'elemento u_ip=1 indica
%      che una serie temporale i e' stata assegnata al cluster p,
%      altrimenti non e' assegnata a p.
% - k: numero di cluster da individuare.
% - m: numero totale di time stamps per ogni serie temporale.
% Output:
% - Z: matrice dei centroidi dei k cluster, rappresentate da
%      vettori Z={Z_1,Z_2,...Z_k) che contengono le serie temporali che
%      meglio rappresentano l’andamento delle serie appartenenti ad un
%      determinato cluster, i centroidi sono quindi serie temporali.

% fissato W0 e U trovo la matrice Z utilizzando la formula (5)
Z = zeros(k, m);
for p = 1:k
    for j = 1:m 
        Z(p, j) = sum(U(:, p).*X(:, j)) / sum(U(:, p));
    end
end
end

function [W] = solve_W(X, U, Z, k, n, m, H, Aeq, beq, lb, ub)
% Author: Elia Mercatanti
% Email: elia.mercatanti@stud.unifi.it
% 
% solve_W(X, U, Z, k, n, m, H, Aeq, beq, lb, ub):
%   Ricava la matrice W minimizzando la funzione obiettivo 1/2*w'*H*w+f'*w,
%   prendendo come costanti le matrici U e Z. Essendo un problema di 
%   programmazione quadratica viene utilizzato il tool matlab 'quadprog'. 
%   Viene inoltre calcolato il vettore f che rappresenta il termine lineare
%   nell'espressione 1/2*w'*H*w+f'*w, necessario per l'uso di quadprog.
% Reference: http://dx.doi.org/10.1016/j.ins.2016.05.040
%
% Input:
% - X: matrice di n serie temporali X={X_1,X_2,...X_n), rappresentate da
%      vettori. Ogni vettore di una serie temporale X_i={x_i1,x_i2,...x_im)
%      e' caratterizzato da m valori che coincidono con m time stamps.
% - U: matrice di appartenenza binaria n x k, dove l'elemento u_ip=1 indica
%      che una serie temporale i e' stata assegnata al cluster p,
%      altrimenti non e' assegnata a p.
% - Z: matrice dei centroidi dei k cluster, rappresentate da vettori 
%      Z={Z_1,Z_2,...Z_k) che contengono le serie temporali che meglio
%      rappresentano l’andamento delle serie appartenenti ad un determinato
%      cluster, i centroidi sono quindi serie temporali.
% - k: numero di cluster da individuare.
% - n: numero di serie temporali.
% - m: numero totale di time stamps per ogni serie temporale.
% - H:  matrice simmetrica che rappresenta la forma quadratica di
%       1/2*w'*H*w+f'*w., utilizzata dal tool 'quadprog'.
% - Aeq: matrice che rappresenta i coefficenti lineari del vincolo
%        Aeq*w=beq.
% - beq: rappresenta il vettore costante nel vincolo Aeq*x = beq.
% - lb: vettore che rappresenta il limite inferiore nel vincolo lb?w?ub.
% - ub: vettore che rappresenta il limite superiore nel vincolo lb?w?ub.
% Output:
% - W: W={W_1,W_2,...W_k) è un set di k vettori che rappresentatno i pesi
%      dei vari time stamps in ogni cluster. Il valore dell'elemento w_pj
%      indica il peso del j-esimo time stamp per il p-esimo cluster.

% ricavo il vettore f per il tool quadprog, rappresenta i coefficenti della
% parte lineare della funzione obiettivo
f = zeros(k*m, 1);
for p = 1:k
    for i = 1:n
        for j = 1:m
            f(j+m*(p-1)) = f(j+m*(p-1)) + U(i, p)*((X(i, j)-Z(p, j))^2);
        end
    end
end

% risolvo il problema di programmazione quadratica e ricostruisco
% correttamente la matrice W
options = optimset('Display', 'off');
W = quadprog(H, f, [], [], Aeq, beq, lb, ub, [], options);
W = vec2mat(W, m);
end

function [x,v] = randfixedsum(n,m,s,a,b)
% [x,v] = randfixedsum(n,m,s,a,b)
%
% This generates an n by m array x, each of whose m columns
% contains n random values lying in the interval [a,b], but
% subject to the condition that their sum be equal to s.  The
% scalar value s must accordingly satisfy n*a <= s <= n*b.  The
% distribution of values is uniform in the sense that it has the
% conditional probability distribution of a uniform distribution
% over the whole n-cube, given that the sum of the x's is s.
%
% The scalar v, if requested, returns with the total
% n-1 dimensional volume (content) of the subset satisfying
% this condition.  Consequently if v, considered as a function
% of s and divided by sqrt(n), is integrated with respect to s
% from s = a to s = b, the result would necessarily be the
% n-dimensional volume of the whole cube, namely (b-a)^n.
%
% This algorithm does no "rejecting" on the sets of x's it
% obtains.  It is designed to generate only those that satisfy all
% the above conditions and to do so with a uniform distribution.
% It accomplishes this by decomposing the space of all possible x
% sets (columns) into n-1 dimensional simplexes.  (Line segments,
% triangles, and tetrahedra, are one-, two-, and three-dimensional
% examples of simplexes, respectively.)  It makes use of three
% different sets of 'rand' variables, one to locate values
% uniformly within each type of simplex, another to randomly
% select representatives of each different type of simplex in
% proportion to their volume, and a third to perform random
% permutations to provide an even distribution of simplex choices
% among like types.  For example, with n equal to 3 and s set at,
% say, 40% of the way from a towards b, there will be 2 different
% types of simplex, in this case triangles, each with its own
% area, and 6 different versions of each from permutations, for
% a total of 12 triangles, and these all fit together to form a
% particular planar non-regular hexagon in 3 dimensions, with v
% returned set equal to the hexagon's area.
%
% Roger Stafford - Jan. 19, 2006

% Check the arguments.
if (m~=round(m)) || (n~=round(n)) || (m<0) || (n<1)
    error('n must be a whole number and m a non-negative integer.')
elseif (s<n*a) || (s>n*b) || (a>=b)
    error('Inequalities n*a <= s <= n*b and a < b must hold.')
end

% Rescale to a unit cube: 0 <= x(i) <= 1
s = (s-n*a)/(b-a);

% Construct the transition probability table, t.
% t(i,j) will be utilized only in the region where j <= i + 1.
k = max(min(floor(s),n-1),0); % Must have 0 <= k <= n-1
s = max(min(s,k+1),k); % Must have k <= s <= k+1
s1 = s - (k:-1:k-n+1); % s1 & s2 will never be negative
s2 = (k+n:-1:k+1) - s;
w = zeros(n,n+1); w(1,2) = realmax; % Scale for full 'double' range
t = zeros(n-1,n);
tiny = 2^(-1074); % The smallest positive matlab 'double' no.
for i = 2:n
    tmp1 = w(i-1,2:i+1).*s1(1:i)/i;
    tmp2 = w(i-1,1:i).*s2(n-i+1:n)/i;
    w(i,2:i+1) = tmp1 + tmp2;
    tmp3 = w(i,2:i+1) + tiny; % In case tmp1 & tmp2 are both 0,
    tmp4 = (s2(n-i+1:n) > s1(1:i)); % then t is 0 on left & 1 on right
    t(i-1,1:i) = (tmp2./tmp3).*tmp4 + (1-tmp1./tmp3).*(~tmp4);
end

% Derive the polytope volume v from the appropriate
% element in the bottom row of w.
v = n^(3/2)*(w(n,k+2)/realmax)*(b-a)^(n-1);

% Now compute the matrix x.
x = zeros(n,m);
if m == 0, return, end % If m is zero, quit with x = []
rt = rand(n-1,m); % For random selection of simplex type
rs = rand(n-1,m); % For random location within a simplex
s = repmat(s,1,m);
j = repmat(k+1,1,m); % For indexing in the t table
sm = zeros(1,m); pr = ones(1,m); % Start with sum zero & product 1
for i = n-1:-1:1  % Work backwards in the t table
    e = (rt(n-i,:)<=t(i,j)); % Use rt to choose a transition
    sx = rs(n-i,:).^(1/i); % Use rs to compute next simplex coord.
    sm = sm + (1-sx).*pr.*s/(i+1); % Update sum
    pr = sx.*pr; % Update product
    x(n-i,:) = sm + pr.*e; % Calculate x using simplex coords.
    s = s - e; j = j - e; % Transition adjustment
end
x(n,:) = sm + pr.*s; % Compute the last x

% Randomly permute the order in the columns of x and rescale.
rp = rand(n,m); % Use rp to carry out a matrix 'randperm'
[~,p] = sort(rp); % The values placed in ig are ignored
x = (b-a)*x(p+repmat((0:n:n*(m-1)),n,1))+a; % Permute & rescale x

return
end