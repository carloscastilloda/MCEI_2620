% =========================================================================
%  Taller de Algebra Lineal Numerica - MCEI
%  PUNTO 1: Inversa de Moore-Penrose (pseudoinversa)
%
%  A es 4x3. Se calcula A^+ y se verifican las 4 propiedades de Penrose:
%     (1) A A^+ A   = A
%     (2) A^+ A A^+ = A^+
%     (3) (A A^+)' = A A^+      (A A^+ simetrica)
%     (4) (A^+ A)' = A^+ A      (A^+ A simetrica)
%
% =========================================================================
clear; close all; format long g;

A = [ 1  0  2;
      2 -1  5;
      0  1 -1;
      1  3 -1];

[m, n] = size(A);
printf('\n PUNTO 1: PSEUDOINVERSA DE MOORE-PENROSE \n');
printf('A (%dx%d) =\n', m, n); disp(A);

% Rango numerico via SVD (no via det, que no aplica a matrices no cuadradas)
[U, S, V] = svd(A);
sv  = diag(S);
tol = max(m, n) * eps(max(sv));          % misma tolerancia que usa pinv()
r   = sum(sv > tol);

printf('\nValores singulares:\n');
for k = 1:numel(sv)
  printf('   sigma_%d = %.12e\n', k, sv(k));
end
printf('Tolerancia de rango  : %.6e\n', tol);
printf('Rango numerico       : %d  (de un maximo de %d)\n', r, min(m,n));
if r < min(m,n)
  printf('>> A es DEFICIENTE DE RANGO: A''A es singular,\n');
  printf('   la formula (A''A)^-1 A'' NO es aplicable. Se usa SVD.\n');
end
printf('Numero de condicion (sigma_1/sigma_r) = %.6e\n', sv(1)/sv(r));

% pseudoinversa via SVD
% A = U S V'   =>   A^+ = V S^+ U'  con S^+ = diag(1/sigma_i) solo para
% los sigma_i > tol (los demas se ponen en cero: truncacion).
Sp = zeros(n, m);
for k = 1:r
  Sp(k,k) = 1 / sv(k);
end
Aplus_svd = V * Sp * U';

% pseudoinversa con la funcion interna de Octave (referencia)
Aplus_oct = pinv(A);

printf('\nA^+ calculada con SVD (manual) =\n'); disp(Aplus_svd);
printf('||A^+_SVD - pinv(A)||_F = %.6e\n', norm(Aplus_svd - Aplus_oct, 'fro'));

Ap = Aplus_svd;   % la que se usa para verificar

% intento por ecuaciones normales
printf('\n Que pasa con la formula (A''A)^{-1} A''? \n');
AtA = A' * A;
printf('det(A''A)  = %.6e   (≈ 0  => singular)\n', det(AtA));
printf('rcond(A''A) = %.6e   (≈ 0  => mal condicionada/singular)\n', rcond(AtA));
warning('off', 'Octave:nearly-singular-matrix');
try
  Ap_ne = inv(AtA) * A';
  printf('||(A''A)^{-1}A'' - A^+||_F = %.6e  <-- basura numerica\n', ...
         norm(Ap_ne - Ap, 'fro'));
catch err
  printf('Fallo: %s\n', err.message);
end
warning('on', 'Octave:nearly-singular-matrix');

% 4 propiedades de Penrose
printf('\n--- VERIFICACION DE LAS 4 PROPIEDADES DE PENROSE ---\n');
e1 = norm(A*Ap*A - A,            'fro');
e2 = norm(Ap*A*Ap - Ap,          'fro');
e3 = norm((A*Ap)' - (A*Ap),      'fro');
e4 = norm((Ap*A)' - (Ap*A),      'fro');

printf('(1) ||A A^+ A - A||_F     = %.6e   %s\n', e1, veredicto(e1));
printf('(2) ||A^+ A A^+ - A^+||_F = %.6e   %s\n', e2, veredicto(e2));
printf('(3) ||(A A^+)'' - A A^+|| = %.6e   %s\n', e3, veredicto(e3));
printf('(4) ||(A^+ A)'' - A^+ A|| = %.6e   %s\n', e4, veredicto(e4));

printf('\nA A^+ (proyector sobre Im(A), traza = rango = %d):\n', r);
disp(A*Ap);  printf('traza = %.12f\n', trace(A*Ap));
printf('\nA^+ A (proyector sobre Im(A''), traza = rango = %d):\n', r);
disp(Ap*A);  printf('traza = %.12f\n', trace(Ap*A));
printf('\nNota: A^+A ~= I porque rank(A)=%d < n=%d. La pseudoinversa\n', r, n);
printf('      NO es inversa por izquierda ni por derecha en este caso.\n');
