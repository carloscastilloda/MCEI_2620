% =========================================================================
%  Taller de Algebra Lineal Numerica - MCEI
%  PUNTO 3: Condicionamiento numerico y estabilidad
%
% =========================================================================
clear; close all; format long g;

printf('\n PUNTO 3: CONDICIONAMIENTO Y ESTABILIDAD \n');

% 3.1 numero de condicion
printf('\n--- 3.1 cond() de la matriz de Hilbert H(i,j)=1/(i+j-1) ---\n');
printf('  n |    cond_2(H)     | digitos perdidos | digitos utiles (de 16)\n');
for n = 2:12
  H  = hilb(n);
  c  = cond(H);
  dp = log10(c);
  printf('%3d | %16.4e | %16.1f | %6.1f\n', n, c, dp, 16-dp);
end

n = 10;                       % caso de estudio
H = hilb(n);
printf('\nMatriz de estudio: H = hilb(%d)\n', n);
printf('cond_2(H)   = %.6e\n', cond(H));
printf('cond_inf(H) = %.6e\n', cond(H, inf));
printf('rcond(H)    = %.6e   (1/rcond ~ cond_1)\n', rcond(H));
sv = svd(H);
printf('Valores singulares (sigma_1 ... sigma_n):\n');
printf('   %.6e\n', sv);
printf('sigma_max/sigma_min = %.6e = cond_2(H)\n', sv(1)/sv(end));

% Comparacion contra una matriz BIEN condicionada (la del punto 2)
Ab = [2 1 0 3 2 1 0 2 1 4; 1 3 2 0 1 4 2 1 0 2; 0 2 4 1 3 0 1 2 4 1;
      3 0 1 5 2 1 3 0 2 1; 2 1 3 2 6 2 1 4 0 3; 1 4 0 1 2 5 2 1 3 0;
      0 2 1 3 1 2 4 0 2 1; 2 1 2 0 4 1 0 5 3 2; 1 0 4 2 0 3 2 3 6 1;
      4 2 1 1 3 0 1 2 1 5];
printf('\nReferencia bien condicionada (A del punto 2): cond_2 = %.4f\n', cond(Ab));

% 3.2 solucion original vs perturbada
x_exact = ones(n,1);
b       = H * x_exact;

delta   = 1e-10;                 % perturbacion pequena y determinista
db      = delta * ones(n,1);
b_pert  = b + db;

x0 = H \ b;                      % sistema original
x1 = H \ b_pert;                 % sistema perturbado en b

rel_b = norm(db)      / norm(b);
rel_x = norm(x1 - x0) / norm(x0);

printf('\n--- 3.2 Perturbacion del lado derecho b ---\n');
printf('Perturbacion absoluta        : ||db||   = %.4e\n', norm(db));
printf('Perturbacion RELATIVA en b   : %.6e\n', rel_b);
printf('Cambio RELATIVO en x         : %.6e\n', rel_x);
printf('Factor de amplificacion real : %.4e\n', rel_x / rel_b);
printf('Cota teorica cond_2(H)       : %.4e   (>= factor real, se cumple)\n', cond(H));

printf('\n  i |   x original   |  x perturbado  |  dif. absoluta\n');
for i = 1:n
  printf('%3d | %14.8f | %14.8f | %14.6e\n', i, x0(i), x1(i), x1(i)-x0(i));
end

% error frente a la solucion exacta conocida (todo unos)
printf('\nError relativo vs x_exacta (sin perturbar) : %.4e\n', ...
       norm(x0 - x_exact)/norm(x_exact));
printf('Error relativo vs x_exacta (perturbado)    : %.4e\n', ...
       norm(x1 - x_exact)/norm(x_exact));
printf('NOTA: incluso SIN perturbar hay error, porque el solo hecho de\n');
printf('      almacenar H en doble precision ya es una perturbacion ~1e-16.\n');

% la misma perturbacion en la matriz BIEN cond.
b_ok  = Ab * x_exact;
x0_ok = Ab \ b_ok;
x1_ok = Ab \ (b_ok + db);
printf('\n--- Contraste con la matriz BIEN condicionada ---\n');
printf('cond=%.2f : pert.rel. b = %.3e -> cambio rel. x = %.3e (amplif. %.2e)\n', ...
       cond(Ab), norm(db)/norm(b_ok), norm(x1_ok-x0_ok)/norm(x0_ok), ...
       (norm(x1_ok-x0_ok)/norm(x0_ok)) / (norm(db)/norm(b_ok)));
printf('La MISMA perturbacion es inofensiva: el condicionamiento es una\n');
printf('propiedad del PROBLEMA, no del algoritmo.\n');

% 3.3 estabilidad de los distintos algoritmos
printf('\n--- 3.3 Estabilidad de los metodos sobre el problema mal condicionado ---\n');

% Gauss/LU
[L, U, P] = lu(H);           x_lu  = U \ (L \ (P*b));
                             x_lu2 = U \ (L \ (P*b_pert));
% QR
[Q, R] = qr(H);              x_qr  = R \ (Q' * b);
                             x_qr2 = R \ (Q' * b_pert);
% SVD (con y sin truncacion)
[Us, Ss, Vs] = svd(H);  s = diag(Ss);
x_svd  = Vs * ((Us' * b)      ./ s);
x_svd2 = Vs * ((Us' * b_pert) ./ s);

% SVD truncada: se descartan los sigma_i < tol (regularizacion)
tolr = 1e-12 * s(1);
k    = sum(s > tolr);
x_tr = Vs(:,1:k) * ((Us(:,1:k)' * b) ./ s(1:k));

printf('%-14s %-14s %-14s %-14s\n', 'Metodo', 'err.rel.x', '||Hx-b||', 'sensib.pert.');
nom = {'LU (Gauss)', 'QR', 'SVD completa'};
Xa  = [x_lu,  x_qr,  x_svd ];
Xb  = [x_lu2, x_qr2, x_svd2];
for j = 1:3
  printf('%-14s %-14.4e %-14.4e %-14.4e\n', nom{j}, ...
         norm(Xa(:,j)-x_exact)/norm(x_exact), ...
         norm(H*Xa(:,j)-b), ...
         norm(Xb(:,j)-Xa(:,j))/norm(Xa(:,j)));
end
printf('%-14s %-14.4e %-14.4e  (rango efectivo usado: %d de %d)\n', ...
       'SVD truncada', norm(x_tr-x_exact)/norm(x_exact), norm(H*x_tr-b), k, n);

printf('\nLECTURA:\n');
printf(' * El RESIDUAL ||Hx-b|| es pequenisimo en los tres metodos: son\n');
printf('   algoritmos ESTABLES (backward stable).\n');
printf(' * El ERROR ||x-x_exacta|| es grande en los tres: el PROBLEMA esta\n');
printf('   mal condicionado. Ningun algoritmo puede arreglar eso.\n');
printf(' * Residual pequeno NO implica solucion correcta cuando cond(A) es alto.\n');

