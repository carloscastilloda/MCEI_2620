% =====================================================================
% Taller integracion y diferenciacion numerica
% MCEI_M - Escuela Colombiana de Ingenieria
% Elaborado por: Carlos Alberto Castillo
% =====================================================================

clear; clc; close all;

function s = ifelse_str(c, a, b)
  if c, s = a; else, s = b; end
end

%%  EJERCICIO 1: integración de una función
f = @(x) exp(-0.4*x) .* (1 + 0.5*sin(3*x));
a = 0;
b = 8;

% --- Referencia analítica ---

al = 0.4; be = 3; L = 8;
I_ref = (1 - exp(-al*L))/al + ...
        0.5*(be - exp(-al*L)*(al*sin(be*L) + be*cos(be*L)))/(al^2 + be^2);
fprintf('Referencia analítica: %.15f\n\n', I_ref);

%  Trapecio compuesto (implementación explícita)
function I = trapecio_comp(f, a, b, n)
  x = linspace(a, b, n+1);
  fx = f(x);
  h = (b - a)/n;
  I = h*(0.5*(fx(1) + fx(end)) + sum(fx(2:end-1)));
end


trap = @(n) trapecio_comp(f, a, b, n);

N = [10 20 50 100 500 1000];
nrep = 2000;                         % repeticiones para medir tiempos
fprintf('%6s %20s %14s %12s %8s\n', 'n', 'Integral trapecio', 'Error abs', 'Tiempo [s]', 'Orden');

res1 = zeros(numel(N), 4);

for k = 1:numel(N)
  n = N(k);
  tic; for r = 1:nrep, I = trap(n); end; t = toc/nrep;
  e = abs(I - I_ref);
  if k > 1, p = log(res1(k-1,3)/e)/log(n/N(k-1)); else, p = NaN; end
  res1(k,:) = [n I e t];
  fprintf('%6d %20.15f %14.3e %12.3e %8.3f\n', n, I, e, t, p);
end

%  Función integrada de Octave: integral (adaptativa) y quadgk
tols = [1e-6 1e-10 1e-14];

fprintf('\n%10s %20s %14s %12s\n', 'RelTol', 'integral()', 'Error abs', 'Tiempo [s]');

for tol = tols
  tic; for r = 1:200, Iq = integral(f, a, b, 'RelTol', tol, 'AbsTol', 1e-15); end; t = toc/200;
  fprintf('%10.0e %20.15f %14.3e %12.3e\n', tol, Iq, abs(Iq - I_ref), t);
end

[Igk, errgk] = quadgk(f, a, b, 'RelTol', 1e-10);
fprintf('quadgk: I = %.15f, error estimado = %.2e, error real = %.2e\n', ...
        Igk, errgk, abs(Igk - I_ref));

% Gráfica: función y área acumulada
xx = linspace(a, b, 1000);
Fcum = cumtrapz(xx, f(xx));

figure(1);
subplot(2,1,1); plot(xx, f(xx), 'b', 'LineWidth', 1.5); grid on;
xlabel('x'); ylabel('f(x)'); title('f(x) = e^{-0.4x}(1+0.5 sin 3x)');
subplot(2,1,2); plot(xx, Fcum, 'r', 'LineWidth', 1.5); grid on;
xlabel('x'); ylabel('F(x) = \int_0^x f'); title('Área acumulada');
print('-dpng', 'ej1_funcion_area.png');

figure(2);
loglog(res1(:,1), res1(:,3), 'o-', 'LineWidth', 1.5); hold on;
loglog(res1(:,1), res1(1,3)*(res1(1,1)./res1(:,1)).^2, 'k--');
grid on; xlabel('n'); ylabel('|error|'); legend('Trapecio', 'O(h^2)');
title('Convergencia del trapecio');
print('-dpng', 'ej1_convergencia.png');

%%  EJERCICIO 2: 50 datos equiespaciados
D = dlmread('/home/carlos/MCEI_2620/Python/jupyter_lab/datos_sensor.csv', ',', 1, 0);
x = D(:,1); y = D(:,2);

% Verificaciones
assert(numel(x) == 50, 'No hay 50 datos');
assert(~any(isnan(D(:))), 'Hay valores faltantes');
h = diff(x);
assert(max(abs(h - h(1))) < 1e-9, 'No es equiespaciado');
h = h(1); n = numel(x) - 1;
fprintf('\nDatos: %d puntos, h = %.4f, %d subintervalos (%s)\n', ...
        numel(x), h, n, ifelse_str(mod(n,2)==0, 'par', 'IMPAR'));

% Trapecio
I_trap = trapz(x, y);

function I = simpson13(y, h)
  % y con número impar de puntos (número par de intervalos)
  n = numel(y) - 1;
  if mod(n, 2) ~= 0, error('Simpson 1/3 requiere n par');
    end
  I = h/3*(y(1) + y(end) + 4*sum(y(2:2:end-1)) + 2*sum(y(3:2:end-2)));
end

function I = simpson38(y, h)
  % exactamente 4 puntos (3 intervalos)
  I = 3*h/8*(y(1) + 3*y(2) + 3*y(3) + y(4));
end



I_simp = simpson13(y(1:47), h) + simpson38(y(47:50), h);

% Alternativa:
I_simp_trap = simpson13(y(1:49), h) + h*(y(49) + y(50))/2;

fprintf('Trapecio                    : %.10f\n', I_trap);
fprintf('Simpson 1/3 (46) + 3/8 (3)  : %.10f\n', I_simp);
fprintf('Simpson 1/3 (48) + trap (1) : %.10f\n', I_simp_trap);
fprintf('Diferencia trapecio-Simpson : %.3e (%.4f %%)\n', ...
        I_trap - I_simp, 100*abs(I_trap - I_simp)/abs(I_simp));

figure(3);
area(x, y, 'FaceColor', [0.7 0.85 1]); hold on;
plot(x, y, 'ko-', 'MarkerFaceColor', 'k', 'MarkerSize', 3);
grid on; xlabel('x'); ylabel('y');
title(sprintf('Datos del sensor - I_{trap} = %.4f,  I_{Simp} = %.4f', I_trap, I_simp));
print('-dpng', 'ej2_datos_area.png');

%%  EXTENSIÓN: diferenciación numérica
% Hacia adelante (O(h)) y centrada (O(h^2)); extremos con fórmulas
% unilaterales de 2º orden para no perder precisión.
dy_fwd = [diff(y)/h; NaN];                    % último punto sin dato
dy_cen = zeros(size(y));
dy_cen(2:end-1) = (y(3:end) - y(1:end-2))/(2*h);
dy_cen(1)   = (-3*y(1) + 4*y(2) - y(3))/(2*h);
dy_cen(end) = ( 3*y(end) - 4*y(end-1) + y(end-2))/(2*h);
dy_grad = gradient(y, h);                     % Octave: extremos O(h)

% Derivada exacta (solo para evaluar el error; no se usa en el cálculo)
dy_ex = 0.35*0.7*cos(0.7*x) - 0.15*2.1*sin(2.1*x) + 0.03;
fprintf('\nError máx. derivada adelante : %.3e\n', max(abs(dy_fwd(1:end-1) - dy_ex(1:end-1))));
fprintf('Error máx. derivada centrada : %.3e\n', max(abs(dy_cen - dy_ex)));

% Sensibilidad al ruido
rand('seed', 1); randn('seed', 1);
sig = 0.01;
yn = y + sig*randn(size(y));
dyn = zeros(size(y));
dyn(2:end-1) = (yn(3:end) - yn(1:end-2))/(2*h);
dyn([1 end]) = NaN;
fprintf('Ruido sigma=%.3f -> cambio integral: %.3e | error máx derivada: %.3e\n', ...
        sig, abs(trapz(x, yn) - I_trap), max(abs(dyn(2:end-1) - dy_ex(2:end-1))));

figure(4);
plot(x, dy_ex, 'k-', x, dy_fwd, 'bs--', x, dy_cen, 'ro-', x, dyn, 'g.-');
grid on; xlabel('x'); ylabel('dy/dx');
legend('Exacta (ref.)', 'Adelante', 'Centrada', 'Centrada c/ruido', 'Location', 'best');
title('Diferenciación numérica');
print('-dpng', 'ext_derivadas.png');
