% =========================================================================
%  Taller de Algebra Lineal Numerica - MCEI
%  PUNTO 4: Inversa de una matriz por 3 vias
%             (a) metodo directo  inv(A)      (LU con pivoteo parcial)
%             (b) descomposicion QR           A^-1 = R^-1 Q'
%             (c) descomposicion SVD          A^-1 = V S^-1 U'
%           y evaluacion de cual es mas estable numericamente.
% =========================================================================
clear; close all; format long g;

printf('\n PUNTO 4: INVERSA POR METODO DIRECTO, QR Y SVD \n');

M1 = [2 1 0 3 2 1 0 2 1 4; 1 3 2 0 1 4 2 1 0 2; 0 2 4 1 3 0 1 2 4 1;
      3 0 1 5 2 1 3 0 2 1; 2 1 3 2 6 2 1 4 0 3; 1 4 0 1 2 5 2 1 3 0;
      0 2 1 3 1 2 4 0 2 1; 2 1 2 0 4 1 0 5 3 2; 1 0 4 2 0 3 2 3 6 1;
      4 2 1 1 3 0 1 2 1 5];
M2 = hilb(10);

casos  = {M1, M2};
nombres = {'A del punto 2 (bien condicionada)', 'Hilbert 10x10 (mal condicionada)'};

for c = 1:2
  A = casos{c};
  n = size(A,1);
  I = eye(n);

  printf('CASO %d: %s\n', c, nombres{c});
  printf('cond_2(A) = %.6e\n', cond(A));

  % (a) inversa directa
  t = tic;
  Ainv_dir = inv(A);
  t_dir = toc(t);

  % (b) inversa via QR
  % A = Q R  =>  A^-1 = R^-1 Q'.  Se resuelve R * Ainv = Q' columna a columna
  % (sustitucion regresiva), nunca formando R^-1 explicitamente.
  t = tic;
  [Q, R] = qr(A);
  Ainv_qr = R \ Q';
  t_qr = toc(t);

  % (c) inversa via SVD
  % A = U S V'  =>  A^-1 = V S^-1 U'
  t = tic;
  [U, S, V] = svd(A);
  s = diag(S);
  Ainv_svd = V * diag(1 ./ s) * U';
  t_svd = toc(t);

  % SVD truncada (pseudoinversa regularizada) con dos tolerancias:
  %   tol_std = la que usa pinv() por defecto -> normalmente NO trunca
  %   tol_agr = tolerancia agresiva           -> regulariza de verdad
  k_std    = sum(s > max(size(A)) * eps(s(1)));
  Ainv_tr  = V(:,1:k_std) * diag(1 ./ s(1:k_std)) * U(:,1:k_std)';
  k_agr    = sum(s > 1e-10 * s(1));
  Ainv_tra = V(:,1:k_agr) * diag(1 ./ s(1:k_agr)) * U(:,1:k_agr)';

  %  comparacion
  printf('\n%-18s %-14s %-14s %-18s %-12s\n', 'Metodo', ...
         '||A*Ai-I||_F', '||Ai*A-I||_F', 'dif.rel vs inv()', 'tiempo[s]');
  met = {'Directo inv()', 'QR', 'SVD'};
  Ai  = {Ainv_dir, Ainv_qr, Ainv_svd};
  tt  = [t_dir, t_qr, t_svd];
  nrm = norm(Ainv_dir, 'fro');
  for j = 1:3
    printf('%-18s %-14.4e %-14.4e %-18.4e %-12.3e\n', met{j}, ...
           norm(A*Ai{j} - I, 'fro'), norm(Ai{j}*A - I, 'fro'), ...
           norm(Ai{j} - Ainv_dir, 'fro')/nrm, tt(j));
  end
  printf('%-18s %-14.4e %-14.4e  (rango efectivo %d/%d)\n', 'SVD trunc (std)', ...
         norm(A*Ainv_tr - I, 'fro'), norm(Ainv_tr*A - I, 'fro'), k_std, n);
  printf('%-18s %-14.4e %-14.4e  (rango efectivo %d/%d)\n', 'SVD trunc (1e-10)', ...
         norm(A*Ainv_tra - I, 'fro'), norm(Ainv_tra*A - I, 'fro'), k_agr, n);

  % cota teorica
  printf('\nCota teorica del error relativo: cond_2(A)*eps = %.4e\n', cond(A)*eps);
  printf('Costo: inv/LU ~ 2n^3,  QR ~ 4n^3/3 + n^3,  SVD ~ 12..22 n^3\n');

  %  prueba de uso real: x = Ainv*b con b = A*ones

  x_ex = ones(n,1);
  bvec = A * x_ex;
  metx = {'Directo inv()', 'QR', 'SVD', 'SVD trunc (std)', 'SVD trunc (1e-10)'};
  Aix  = {Ainv_dir, Ainv_qr, Ainv_svd, Ainv_tr, Ainv_tra};
  printf('\nError al RESOLVER A x = b con cada inversa (b = A*ones):\n');
  for j = 1:5
    printf('   %-18s err.rel(x) = %.4e\n', metx{j}, ...
           norm(Aix{j}*bvec - x_ex)/norm(x_ex));
  end

  if c == 2
    printf('\nValores singulares de Hilbert(10):\n'); printf('   %.6e\n', s);
    printf('Los ultimos sigma estan al nivel del ruido de redondeo: al\n');
    printf('invertirlos (1/sigma) el ruido se amplifica ~1e13 veces.\n');
  end
end

%  mensaje final
printf('\n--- CONCLUSION DEL PUNTO 4 ---\n');
printf('1) Si la matriz esta BIEN condicionada los tres metodos dan\n');
printf('   practicamente el mismo resultado (~1e-15); gana inv() por costo.\n');
printf('2) Si esta MAL condicionada NINGUNO se salva: los tres errores\n');
printf('   quedan del orden de cond(A)*eps, que es la cota teorica. La\n');
printf('   perdida de precision es del PROBLEMA, no del algoritmo.\n');
printf('   La ventaja real de la SVD no es un error menor, sino que es el\n');
printf('   unico metodo que DIAGNOSTICA (muestra los sigma pequenos) y\n');
printf('   permite REGULARIZAR truncandolos (ultima fila de la tabla).\n');
printf('3) REGLA DE ORO: para RESOLVER A x = b nunca se calcula A^-1.\n');
printf('   x = A\\b es mas rapido (O(n^3) vs 2n^3 + n^2) y mas preciso.\n');
printf('   La inversa explicita solo se justifica si se necesita la matriz\n');
printf('   en si (p.ej. la matriz de covarianza (X''X)^-1 en estadistica).\n');

