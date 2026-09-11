% =========================================================================
%  Taller de Algebra Lineal Numerica - MCEI
%  PROBLEMA PRACTICO: convertidor DC-DC en regimen estacionario
%
%     Modelo:  P = b0 + b1*V + b2*I + b3*T + e
%     Forma matricial:  y = X*beta,   X = [1 V I T]  (n x 4)
%
% =========================================================================
clear; close all; format long g;

archivo = '/home/carlos/MCEI_2620/C_C++/matrices/datos_convertidor_realista.csv';

printf('\n PROBLEMA PRACTICO: CONVERTIDOR DC-DC \n');

% carga de datos
% El CSV trae encabezado: Voltaje_V,Corriente_A,Temperatura_C,Potencia_W
if ~exist(archivo, 'file')
  error(['No se encontro el archivo de datos:\n  %s\n' ...
         'Ejecute este script desde la carpeta octave/ del taller.'], archivo);
end
try
  D = dlmread(archivo, ',', 1, 0);        % se salta 1 fila de encabezado
catch
  % respaldo por si dlmread no esta disponible en su version de Octave
  fid = fopen(archivo, 'r');
  fgetl(fid);                             % descarta el encabezado
  C   = textscan(fid, '%f%f%f%f', 'Delimiter', ',');
  fclose(fid);
  D = cell2mat(C);
end

V = D(:,1);      % voltaje de entrada [V]
Ic = D(:,2);     % corriente de carga [A]   (Ic para no chocar con eye/I)
T = D(:,3);      % temperatura del disipador [C]
y = D(:,4);      % potencia de salida [W]
n = numel(y);

printf('Muestras leidas: n = %d\n', n);
printf('\n%-16s %10s %10s %10s %10s\n', 'Variable','media','desv.est','min','max');
etq = {'Voltaje [V]','Corriente [A]','Temp [C]','Potencia [W]'};
for k = 1:4
  printf('%-16s %10.4f %10.4f %10.4f %10.4f\n', etq{k}, ...
         mean(D(:,k)), std(D(:,k)), min(D(:,k)), max(D(:,k)));
end

% matriz de diseno
X = [ones(n,1), V, Ic, T];            % n x 4
p = size(X,2);

printf('\nMatriz de diseno X: %d x %d\n', n, p);
printf('Primeras 3 filas:\n'); disp(X(1:3,:));

% (a) ECUACIONES NORMALES
XtX   = X' * X;
Xty   = X' * y;
t = tic;
b_ne  = XtX \ Xty;            % se resuelve el sistema, no se invierte
t_ne  = toc(t);
b_ne_inv = inv(XtX) * Xty;    % version "literal" con inversa explicita


% (b) QR  (A = QR delgada; R*beta = Q'*y)
t = tic;
[Q, R] = qr(X, 0);            % QR reducida: Q es n x p, R es p x p
b_qr   = R \ (Q' * y);
t_qr   = toc(t);

% (c) SVD / pseudoinversa
t = tic;
[U, S, Vs] = svd(X, 0);       % SVD economica
s          = diag(S);
b_svd      = Vs * ((U' * y) ./ s);
t_svd      = toc(t);


% (d) BACKSLASH (referencia de Octave: QR con pivoteo de columnas)
t = tic;  b_bs = X \ y;  t_bs = toc(t);

%  ENTREGABLE 1: coeficientes
printf('\n===== 1. COEFICIENTES ESTIMADOS =====\n');
printf('%-12s %14s %14s %14s %14s\n', 'coef', 'Ec.normales', 'QR', 'SVD', 'backslash');
nomb = {'b0 (interc)','b1 (V)','b2 (I)','b3 (T)'};
for k = 1:p
  printf('%-12s %14.8f %14.8f %14.8f %14.8f\n', nomb{k}, ...
         b_ne(k), b_qr(k), b_svd(k), b_bs(k));
end

printf('\nDiferencias respecto a QR (norma 2):\n');
printf('  Ec. normales (\\)   : %.4e\n', norm(b_ne     - b_qr));
printf('  Ec. normales inv() : %.4e\n', norm(b_ne_inv - b_qr));
printf('  SVD                : %.4e\n', norm(b_svd    - b_qr));
printf('  backslash          : %.4e\n', norm(b_bs     - b_qr));
printf('Tiempos [s]: normales %.3e | QR %.3e | SVD %.3e | \\ %.3e\n', ...
       t_ne, t_qr, t_svd, t_bs);

beta = b_qr;                  % se adopta la solucion por QR
printf('\n>>> MODELO AJUSTADO:\n');
printf('    P = %.4f %+.4f*V %+.4f*I %+.4f*T\n', ...
       beta(1), beta(2), beta(3), beta(4));

%  ENTREGABLE 3: residual y errores
r    = y - X*beta;
SSE  = r' * r;
MSE  = SSE / n;
RMSE = sqrt(MSE);
MAE  = mean(abs(r));
SST  = sum((y - mean(y)).^2);
R2   = 1 - SSE/SST;
R2a  = 1 - (1-R2)*(n-1)/(n-p);

printf('\n===== 3. RESIDUAL Y ERROR CUADRATICO MEDIO =====\n');
printf('||r||_2 (norma del residual) = %.6f  W\n', norm(r));
printf('SSE  (suma cuadrados error)  = %.6f  W^2\n', SSE);
printf('MSE  = SSE/n                 = %.6f  W^2\n', MSE);
printf('RMSE = sqrt(MSE)             = %.6f  W\n', RMSE);
printf('MAE  (error absoluto medio)  = %.6f  W\n', MAE);
printf('R^2                          = %.8f\n', R2);
printf('R^2 ajustado                 = %.8f\n', R2a);
printf('media(r) = %.4e  (debe ser ~0)   max|r| = %.4f W\n', mean(r), max(abs(r)));
printf('RMSE / media(P) = %.2f %%\n', 100*RMSE/mean(y));
printf('Ortogonalidad ||X''r||_inf    = %.4e  (condicion de MC, debe ser ~0)\n', ...
       norm(X'*r, inf));

%  ENTREGABLE 2: influencia relativa de cada variable
sd_x = [std(V), std(Ic), std(T)];
sd_y = std(y);
beta_std = beta(2:4)' .* sd_x / sd_y;       % coeficientes estandarizados (beta)
contrib  = 100 * abs(beta_std) / sum(abs(beta_std));
rango    = beta(2:4)' .* [max(V)-min(V), max(Ic)-min(Ic), max(T)-min(T)];

% errores estandar e inferencia
sigma2 = SSE / (n - p);
C      = sigma2 * inv(XtX);                 % matriz de covarianza de beta
se     = sqrt(diag(C));
tval   = beta ./ se;

printf('\n===== 2. INFLUENCIA RELATIVA DE CADA VARIABLE =====\n');
printf('%-10s %12s %12s %12s %12s %12s\n', 'var', 'beta', 'err.est', 't', ...
       'beta_std', 'influencia%');
for k = 1:3
  printf('%-10s %12.6f %12.6f %12.4f %12.6f %12.2f\n', ...
         nomb{k+1}, beta(k+1), se(k+1), tval(k+1), beta_std(k), contrib(k));
end
printf('%-10s %12.6f %12.6f %12.4f %12s %12s\n', nomb{1}, beta(1), se(1), tval(1), '-', '-');

printf('\nEfecto sobre P al recorrer TODO el rango de cada variable:\n');
printf('  V: %+8.2f W   (rango %.2f V)\n', rango(1), max(V)-min(V));
printf('  I: %+8.2f W   (rango %.2f A)\n', rango(2), max(Ic)-min(Ic));
printf('  T: %+8.2f W   (rango %.2f C)\n', rango(3), max(T)-min(T));

printf('\nMatriz de correlacion [V I T P]:\n');
disp(corr(D));

%  ENTREGABLE 5: numero de condicion de X
printf('\n===== 5. CONDICIONAMIENTO DE X =====\n');
printf('Valores singulares de X: '); printf('%.6f  ', s); printf('\n');
printf('cond_2(X)     = %.6f\n', cond(X));
printf('cond_2(X''X)   = %.6f   (= cond(X)^2)\n', cond(XtX));
printf('cond(X)^2     = %.6f   <-- comprobacion\n', cond(X)^2);
printf('\nEsto es EXACTAMENTE por que QR es preferible a las ecuaciones\n');
printf('normales: formar X''X ELEVA AL CUADRADO el numero de condicion.\n');
printf('Digitos perdidos con QR            : %.2f\n', log10(cond(X)));
printf('Digitos perdidos con ec. normales  : %.2f\n', log10(cond(XtX)));

% centrado de las columnas: mejora enorme del condicionamiento
Xc = [ones(n,1), V-mean(V), Ic-mean(Ic), T-mean(T)];
printf('\ncond_2(X) centrando las columnas   = %.4f  (era %.2f)\n', ...
       cond(Xc), cond(X));
printf('El mal condicionamiento venia de la ESCALA dispar de las columnas\n');
printf('(la de unos vs T ~ 50), no de colinealidad real entre V, I y T.\n');

%  diagnostico fisico extra
printf('\n===== DIAGNOSTICO ADICIONAL: el modelo fisico =====\n');
Pvi  = V .* Ic;
Xvi  = [ones(n,1), Pvi];
bvi  = Xvi \ y;
rvi  = y - Xvi*bvi;
R2vi = 1 - (rvi'*rvi)/SST;
printf('corr(P, V*I) = %.6f\n', corr(y, Pvi));
printf('Ajuste P = a + b*(V*I):  a = %.4f,  b = %.6f\n', bvi(1), bvi(2));
printf('   RMSE = %.4f W   R^2 = %.6f  <-- mucho mejor que el modelo lineal\n', ...
       sqrt((rvi'*rvi)/n), R2vi);
printf('Interpretacion: la potencia real es P = eta*V*I (producto, no suma).\n');
printf('El modelo lineal es solo una LINEALIZACION alrededor del punto de\n');
printf('operacion medio; b ~ %.3f es la eficiencia media del convertidor.\n', bvi(2));

%  efecto del tamano de muestra
printf('\n===== EFECTO DEL NUMERO DE MUESTRAS =====\n');
printf('%6s %10s %10s %10s %10s %10s %10s\n', 'n', 'b0','b1(V)','b2(I)','b3(T)','RMSE','cond(X)');
for m = [10 20 50 100]
  Xm = X(1:m,:);  ym = y(1:m);
  bm = Xm \ ym;   rm = ym - Xm*bm;
  printf('%6d %10.4f %10.4f %10.4f %10.4f %10.4f %10.2f\n', ...
         m, bm(1), bm(2), bm(3), bm(4), sqrt((rm'*rm)/m), cond(Xm));
end
printf('Los coeficientes se estabilizan y los errores estandar caen como\n');
printf('1/sqrt(n); el RMSE tiende al ruido real y cond(X) se estabiliza.\n');

%  graficas
try
  figure();
  subplot(2,2,1); plot(y, X*beta, 'o', [min(y) max(y)], [min(y) max(y)], 'r-');
    xlabel('P medida [W]'); ylabel('P predicha [W]'); title('Ajuste'); grid on;
  subplot(2,2,2); plot(X*beta, r, 'o'); hold on; plot(xlim, [0 0], 'r-');
    xlabel('P predicha [W]'); ylabel('residual [W]'); title('Residuales'); grid on;
  subplot(2,2,3); hist(r, 15); xlabel('residual [W]'); title('Histograma residuales');
  subplot(2,2,4); plot(Pvi, y, 'o'); xlabel('V*I [W]'); ylabel('P [W]');
    title('P vs V*I'); grid on;
  print('-dpng', '-r120', 'ajuste_convertidor.png');
  printf('\n[grafica guardada en ajuste_convertidor.png]\n');
catch
  printf('\n[graficas omitidas: sin backend grafico disponible]\n');
end


