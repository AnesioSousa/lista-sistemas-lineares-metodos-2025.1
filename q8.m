clc
clear all
close all

function [x, k, Erx] = metodo_gauss_seidel(A, b, tol, N, x0)
  % Inicializações
  [n, ~] = size(A);
  x = x0;
  Erx = inf;
  k = 0;

  % Iterações
  while k < N && max(Erx) > tol
    x_old = x;
    for i = 1:n
      sum = 0;
      for j = 1:n
        if j ~= i
          sum += A(i,j) * x(j);
        end
      end
      x(i) = (b(i) - sum) / A(i,i);
    end

    % Calcula Erro Relativo
    Erx = abs((x - x_old) ./ x) * 100;
    k++;
  endwhile

  % Verificação de Convergência
  if max(Erx) > tol
    warning('O método não convergiu');
  endif
endfunction

function beta = criterio_sassenfeld(A)
  n = size(A, 1);
  beta = zeros(n, 1);

  for i = 1:n
    sum = 0;
    for j = 1:i-1
      sum += abs(A(i,j)) * beta(j);
    end
    for j = i+1:n
      sum += abs(A(i,j));
    end
    beta(i) = sum / abs(A(i,i));
  end
endfunction

A = [0.8 -0.4, 0; -0.4, 0.8 -0.4; 0, -0.4 0.8];
b = [41; 25; 105];
x0 = [0; 0; 0];
tol = 0.05;
N = 100;


% Chama o método de Gauss-Seidel
[x, k, Erx] = metodo_gauss_seidel(A, b, tol, N, x0)

% Verifica o critério de Sassenfeld
beta = criterio_sassenfeld(A);
disp('Critério de Sassenfeld:');
disp(beta);
disp('Satisfaz o critério de Sassenfeld?');
disp(all(beta < 1));

% Exibe os resultados
disp('Solução:');
disp(x);
disp('Número de iterações:');
disp(k);
disp('Erro relativo:');
disp(Erx);
