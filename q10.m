clc
clear all
close all

function [x, k, Erx, erros_por_iteracao] = metodo_sor(A, b, tol, N, x0, omega)
  [n, ~] = size(A);
  x = x0;
  Erx = inf;
  k = 0;
  erros_por_iteracao = [];

  while k < N && max(Erx) > tol
    x_old = x;
    for i = 1:n
      soma = 0;
      for j = 1:n
        if j ~= i
          soma += A(i,j) * x(j);
        end
      end
      x_new = (b(i) - soma) / A(i,i);
      x(i) = omega * x_new + (1 - omega) * x(i);
    end

    % Calcula erro relativo e armazena
    Erx = abs((x - x_old) ./ x) * 100;
    erros_por_iteracao = [erros_por_iteracao; max(Erx)];  % guarda o erro máximo da iteração
    k++;
  endwhile

  if max(Erx) > tol
    warning('O método não convergiu');
  endif
endfunction

% Dados do sistema
A = [-8,  1, -2;
      2, -6, -1;
     -3, -1,  7];

b = [-20; -38; -34];
x0 = [0; 0; 0];
tol = 0.05;
omega = 0.95;  % pode ajustar aqui
N = 100;

%omegas = 0.1:0.05:1.95;
%melhor_omega = 0;
%menor_k = Inf;

%for w = omegas
%  [~, k, ~] = metodo_sor(A, b, tol, N, x0, w);
%  if k < menor_k
%    menor_k = k;
%    melhor_omega = w;
%  end
%end

%fprintf('Melhor omega: %.2f\n', melhor_omega);
%fprintf('Menor número de iterações: %d\n', menor_k);

[x, k, Erx, erros] = metodo_sor(A, b, tol, N, x0, omega);

% Resultados
disp('Solução:');
disp(x);
disp('Número de iterações:');
disp(k);

% Gráfico de convergência
figure;
plot(1:length(erros), erros, '-o', 'LineWidth', 2);
xlabel('Iterações');
ylabel('Erro Relativo (%)');
title(['Convergência do Método SOR (\omega = ' num2str(omega) ')']);
grid on;

