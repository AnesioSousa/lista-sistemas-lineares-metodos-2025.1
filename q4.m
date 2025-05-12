clc
clear all
close all

function [L, U, P] = lu_decomposition_pivot(A)
  [n, m] = size(A);
  L = eye(n);
  U = A;
  P = eye(n);

  for j = 1:n-1
    % Mostrar ANTES do pivoteamento
    fprintf('---------------------------------\n');
    fprintf('Passo %d: ANTES do pivoteamento\n', j);
    fprintf('P =\n');
    disp(sprintf('%.4f %.4f %.4f\n', P'));
    fprintf('L =\n');
    disp(sprintf('%.4f %.4f %.4f\n', L'));
    fprintf('U =\n');
    disp(sprintf('%.4f %.4f %.4f\n', U'));

    % Pivoteamento parcial
    [~, pivot_row] = max(abs(U(j:n,j)));
    pivot_row = pivot_row + j - 1;

    if pivot_row ~= j
      % Trocar linhas de U
      temp = U(j,:);
      U(j,:) = U(pivot_row,:);
      U(pivot_row,:) = temp;

      % Trocar linhas de P
      temp = P(j,:);
      P(j,:) = P(pivot_row,:);
      P(pivot_row,:) = temp;

      % Trocar linhas de L (até a coluna j-1)
      if j >= 2
        temp = L(j,1:j-1);
        L(j,1:j-1) = L(pivot_row,1:j-1);
        L(pivot_row,1:j-1) = temp;
      end

      % Mostrar DEPOIS do pivoteamento (somente se houve troca)
      fprintf('Passo %d: DEPOIS do pivoteamento\n', j);
      fprintf('P =\n');
      disp(sprintf('%.4f %.4f %.4f\n', P'));
      fprintf('L =\n');
      disp(sprintf('%.4f %.4f %.4f\n', L'));
      fprintf('U =\n');
      disp(sprintf('%.4f %.4f %.4f\n', U'));
    end

    % Eliminação
    for i = j+1:n
      L(i,j) = U(i,j) / U(j,j);
      U(i,:) = U(i,:) - L(i,j) * U(j,:);
    end

    % Mostrar depois da eliminação
    fprintf('Passo %d: Após ELIMINAÇÃO (Atualizando L e U)\n', j);
    fprintf('L =\n');
    disp(sprintf('%.4f %.4f %.4f\n', L'));
    fprintf('U =\n');
    disp(sprintf('%.4f %.4f %.4f\n', U'));
  end
end

function x = forward_substitution(L, b, show_steps)
  if nargin < 3
    show_steps = true;
  end

  [n, ~] = size(L);
  x = zeros(n, 1);

  for i = 1:n
    x(i) = (b(i) - L(i,1:i-1) * x(1:i-1)) / L(i,i);
    if show_steps
      fprintf('Passo %d: Substituição direta (y)\n', i);
      fprintf('y = [%.4f %.4f %.4f]\n', x(1), x(2), x(3));
    end
  end
end

function x = backward_substitution(U, b, show_steps)
  if nargin < 3
    show_steps = true;
  end

  [n, ~] = size(U);
  x = zeros(n, 1);

  for i = n:-1:1
    x(i) = (b(i) - U(i,i+1:n) * x(i+1:n)) / U(i,i);
    if show_steps
      fprintf('Passo %d: Substituição retroativa (x)\n', i);
      fprintf('x = [%.4f %.4f %.4f]\n', x(1), x(2), x(3));
    end
  end
end

% ----------------- Código Principal ------------------

A = [7  2 -3;
     2  5 -3;
     1  -1 -6];

b = [-12; -20; -26];
b2 = [12; 18; -6];

fprintf('\n\n');
disp('==================================');
disp('Decomposição LU com pivoteamento');
disp('==================================');

% Decomposição LU com pivoteamento
[L, U, P] = lu_decomposition_pivot(A);

fprintf('\n\n');
disp('============================================');
disp('Resolvendo sistema para b = [-12; -20; -26]');
disp('============================================');

% Resolver o sistema para b = [2; 0; 1]
b_mod = P * b;
disp('Vetor Pb:');
fprintf('Pb = [%.4f; %.4f; %.4f]\n', b_mod(1), b_mod(2), b_mod(3));

% Substituição direta
y = forward_substitution(L, b_mod, true);

% Substituição retroativa
x = backward_substitution(U, y, true);

fprintf('\n\n');
disp('================================================');
disp('Agora resolvendo sistema para b2 = [12; 18; -6]');
disp('================================================');

% Resolver o sistema para b2 = [3; 1; 2]
b2_mod = P * b2;
disp('Vetor Pb2:');
fprintf('Pb2 = [%.4f; %.4f; %.4f]\n', b2_mod(1), b2_mod(2), b2_mod(3));

% Substituição direta
y2 = forward_substitution(L, b2_mod, true);

% Substituição retroativa
x2 = backward_substitution(U, y2, true);

A_inv = inv(A);  % <- Adicione essa linha

fprintf('\nNorma de A: %.2f\n', norm(A,inf));
fprintf('Norma da inversa de A: %.2f\n', norm(A_inv,inf));
fprintf('Número Condição:\n');
condInf = norm(A,inf)*norm(A_inv,inf);
fprintf('CondInf: %.2f\n', condInf);

fprintf('\nSolução final para b = [-12; -20; -26]:\n');
fprintf('x = [%.4f; %.4f; %.4f]\n', x(1), x(2), x(3));

fprintf('\nSolução final para b2 = [12; 18; -6]:\n');
fprintf('x = [%.4f; %.4f; %.4f]\n', x2(1), x2(2), x2(3));
