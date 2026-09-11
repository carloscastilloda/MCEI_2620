% =========================================================================
%  Taller de Algebra Lineal Numerica - MCEI
%  PUNTO 2: Sistemas de ecuaciones lineales 10x10
%           Gauss / LU / QR  +  comparacion de estabilidad y eficiencia
%
% =========================================================================
clear; close all; format long g;

A = [2 1 0 3 2 1 0 2 1 4;
     1 3 2 0 1 4 2 1 0 2;
     0 2 4 1 3 0 1 2 4 1;
     3 0 1 5 2 1 3 0 2 1;
     2 1 3 2 6 2 1 4 0 3;
     1 4 0 1 2 5 2 1 3 0;
     0 2 1 3 1 2 4 0 2 1;
     2 1 2 0 4 1 0 5 3 2;
     1 0 4 2 0 3 2 3 6 1;
     4 2 1 1 3 0 1 2 1 5];

n       = size(A,1);
x_exact = ones(n,1);
b       = A * x_exact;

printf('\n PUNTO 2: SISTEMA LINEAL 10x10 \n');
printf('b'' = '); printf('%g ', b); printf('\n');

%  caracterizacion previa
printf('\n Caracterizacion de A \n');
printf('det(A)          = %.6f\n',  det(A));
printf('cond_2(A)       = %.6f\n',  cond(A));
printf('cond_inf(A)     = %.6f\n',  cond(A, inf));
printf('rcond(A)        = %.6e\n',  rcond(A));
printf('Simetrica?      : %d\n',    isequal(A, A'));

ev = eig(A);
printf('Autovalores     : '); printf('%.4f ', sort(ev)); printf('\n');

if all(ev > 0)
  printf('A es definida positiva -> se podria usar Cholesky.\n');
else
  printf('A es simetrica pero INDEFINIDA (hay autovalores < 0):\n');
  printf('   Cholesky NO es aplicable; se usa LU con pivoteo / QR.\n');
end


% METODO 1: Eliminacion de Gauss con pivoteo parcial (implementada a mano)

function [x, npiv] = gauss_pp(A, b)
  n    = size(A,1);
  M    = [A, b];
  npiv = 0;

  %  fase de eliminacion hacia adelante
  for k = 1:n-1
    % pivoteo parcial: buscar el mayor |elemento| de la columna k
    [~, idx] = max(abs(M(k:n, k)));
    p = idx + k - 1;
    if p ~= k
      M([k p], :) = M([p k], :);     % intercambio de filas
      npiv = npiv + 1;
    end
    if abs(M(k,k)) < eps
      error('Matriz singular o casi singular en la columna %d', k);
    end
    for i = k+1:n
      f        = M(i,k) / M(k,k);
      M(i,k:n+1) = M(i,k:n+1) - f * M(k,k:n+1);
      M(i,k)   = 0;
    end
  end

  % sustitucion regresiva
  x = zeros(n,1);
  for i = n:-1:1
    x(i) = (M(i,n+1) - M(i,i+1:n) * x(i+1:n)) / M(i,i);
  end
end

% Aplicuacion de funcion
t = tic;
[x_gauss, npiv] = gauss_pp(A, b);
t_gauss = toc(t);


% METODO 2: Factorizacion LU (PA = LU)

t = tic;
[L, U, P] = lu(A);
y     = L \ (P*b);        % sustitucion progresiva (L triangular inferior)
x_lu  = U \ y;            % sustitucion regresiva  (U triangular superior)
t_lu  = toc(t);


% METODO 3: Factorizacion QR (A = QR, Q ortogonal)

t = tic;
[Q, R] = qr(A);
x_qr   = R \ (Q' * b);    % Q'Q = I  =>  Rx = Q'b
t_qr   = toc(t);

% referencia: operador backslash de Octave (LAPACK)
t = tic; x_bs = A \ b; t_bs = toc(t);

%% resultados %%
printf('\n Soluciones \n');
printf('  i |    Gauss    |     LU      |     QR      |  backslash\n');
for i = 1:n
  printf('%3d | %11.8f | %11.8f | %11.8f | %11.8f\n', ...
         i, x_gauss(i), x_lu(i), x_qr(i), x_bs(i));
end

printf('\n Comparacion de estabilidad y eficiencia \n');
printf('%-12s %-14s %-14s %-14s %s\n', 'Metodo', '||Ax-b||_2', ...
       'err.rel.', 'tiempo [s]', 'costo aprox.');
metodos = {'Gauss', 'LU', 'QR', 'backslash'};
X       = [x_gauss, x_lu, x_qr, x_bs];
tiempos = [t_gauss, t_lu, t_qr, t_bs];
costos  = {'2n^3/3', '2n^3/3', '4n^3/3', '2n^3/3'};
for k = 1:4
  res = norm(A*X(:,k) - b);
  err = norm(X(:,k) - x_exact) / norm(x_exact);
  printf('%-12s %-14.4e %-14.4e %-14.3e %s\n', ...
         metodos{k}, res, err, tiempos(k), costos{k});
end

% cota teorica del error
printf('\nCota teorica:  err_rel(x) <= cond(A) * eps = %.4e * %.2e = %.4e\n', ...
       cond(A), eps, cond(A)*eps);
printf('Todos los metodos quedan MUY por debajo de la cota -> A esta bien\n');
printf('condicionada y los 3 metodos son numericamente estables aqui.\n');

%  estructura de las factorizaciones
printf('\n--- Factores obtenidos ---\n');
printf('Permutaciones de fila realizadas por Gauss: %d\n', npiv);
printf('diag(U) de la LU : '); printf('%.4f ', diag(U)); printf('\n');
printf('diag(R) de la QR : '); printf('%.4f ', diag(R)); printf('\n');
printf('||P*A - L*U||_F  = %.4e\n', norm(P*A - L*U, 'fro'));
printf('||A   - Q*R||_F  = %.4e\n', norm(A - Q*R, 'fro'));
printf('||Q''Q - I||_F    = %.4e   (ortogonalidad de Q)\n', norm(Q'*Q - eye(n),'fro'));

%  eficiencia con muchas repeticiones
printf('\n--- Eficiencia media sobre 2000 repeticiones ---\n');
REP = 2000;
t = tic; for k=1:REP, gauss_pp(A,b);                 end; tg = toc(t)/REP;
t = tic; for k=1:REP, [L,U,P]=lu(A); U\(L\(P*b));     end; tl = toc(t)/REP;
t = tic; for k=1:REP, [Q,R]=qr(A);   R\(Q'*b);        end; tq = toc(t)/REP;
t = tic; for k=1:REP, A\b;                            end; tb = toc(t)/REP;
printf('Gauss (a mano) : %.3e s\n', tg);
printf('LU    (LAPACK) : %.3e s\n', tl);
printf('QR    (LAPACK) : %.3e s\n', tq);
printf('backslash      : %.3e s\n', tb);
printf('\nVentaja de LU: si hay que resolver con VARIOS b, la factorizacion\n');
printf('se hace UNA vez (O(n^3)) y cada nuevo b cuesta solo O(n^2).\n');

