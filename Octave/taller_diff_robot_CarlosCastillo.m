% Parte 1 - GNU Octave: diferencias centrales explícitas
clear; close all; clc;

% Genera trayectoria_robot.csv con columnas t,x,y
% 51 muestras equiespaciadas: t_i = 0.2 i, i = 0,...,50
i = (0:50)';
t = 0.2 * i;

x = 0.08 * t.^2 + 0.40 * sin(0.45 * t);
y = 0.50 * t + 0.30 * (1 - cos(0.45 * t));

fid = fopen('trayectoria_robot.csv', 'w');
fprintf(fid, 't,x,y\n');
fprintf(fid, '%.10f,%.10f,%.10f\n', [t x y]');
fclose(fid);

%% Carga de datos generados
data = dlmread('trayectoria_robot.csv', ',', 1, 0);
t = data(:,1); x = data(:,2); y = data(:,3);
N = numel(t);

tic;
% --- Derivadas de posición: centradas en el interior, unilaterales O(h) en extremos
vx = zeros(N,1); vy = zeros(N,1);
vx(2:N-1) = (x(3:N) - x(1:N-2)) ./ (t(3:N) - t(1:N-2));
vy(2:N-1) = (y(3:N) - y(1:N-2)) ./ (t(3:N) - t(1:N-2));
vx(1) = (x(2)-x(1))/(t(2)-t(1));   vx(N) = (x(N)-x(N-1))/(t(N)-t(N-1));
vy(1) = (y(2)-y(1))/(t(2)-t(1));   vy(N) = (y(N)-y(N-1))/(t(N)-t(N-1));

% --- Velocidad lineal y orientación (con desenvolvimiento)
v     = sqrt(vx.^2 + vy.^2);
theta = unwrap(atan2(vy, vx));

% --- Velocidad angular
omega = zeros(N,1);
omega(2:N-1) = (theta(3:N) - theta(1:N-2)) ./ (t(3:N) - t(1:N-2));
omega(1) = (theta(2)-theta(1))/(t(2)-t(1));
omega(N) = (theta(N)-theta(N-1))/(t(N)-t(N-1));
t_exec = toc;

% --- Solución analítica para validar
xd = 0.16*t + 0.18*cos(0.45*t);          yd = 0.50 + 0.135*sin(0.45*t);
xdd = 0.16 - 0.081*sin(0.45*t);          ydd = 0.06075*cos(0.45*t);
v_ex = sqrt(xd.^2 + yd.^2);
w_ex = (xd.*ydd - yd.*xdd) ./ (xd.^2 + yd.^2);

int = 2:N-1;
printf('Octave | tiempo: %.3e s\n', t_exec);
printf('Error max interior  v: %.3e   omega: %.3e\n', max(abs(v(int)-v_ex(int))), max(abs(omega(int)-w_ex(int))));
printf('Error max global    v: %.3e   omega: %.3e\n', max(abs(v-v_ex)), max(abs(omega-w_ex)));

dlmwrite('resultados_octave.csv', [t vx vy v theta omega], 'precision', '%.10f');

% --- Gráficas
f = figure('position', [100 100 1000 750]);   % visible en pantalla

subplot(2,2,1); plot(x, y, 'o-', 'markersize', 3); axis equal; grid on;
xlabel('x (m)'); ylabel('y (m)'); title('Trayectoria');

subplot(2,2,2); plot(t, v, 'o', 'markersize', 3, t, v_ex, '-'); grid on;
xlabel('t (s)'); ylabel('v (m/s)'); title('Velocidad lineal');
legend('numérica', 'analítica', 'location', 'northwest');

subplot(2,2,3); plot(t, theta, 'o-', 'markersize', 3); grid on;
xlabel('t (s)'); ylabel('\theta (rad)'); title('Orientación (unwrap)');

subplot(2,2,4); plot(t, omega, 'o', 'markersize', 3, t, w_ex, '-'); grid on;
xlabel('t (s)'); ylabel('\omega (rad/s)'); title('Velocidad angular');
legend('numérica', 'analítica', 'location', 'southeast');

drawnow;
print(f, 'graficas_robot.png', '-dpng', '-r110');
printf('Gráfica guardada en: %s\n', fullfile(pwd, 'graficas_octave.png'));
