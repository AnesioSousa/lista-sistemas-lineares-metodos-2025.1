clc
clear all
close all

function [A_inv, b_solution] = gauss_jordan(A, b)
  [n, m] = size(A);
  if n ~= m
    error('A matrix deve ser quadrada');
  end

  M = [A eye(n) b];
  disp('Matriz Inicial Aumentada [A|I|b]:');
  disp(M);

  % Gauss-Jordan Elimination
  for col = 1:n
    fprintf('\n=== Coluna %d ===\n', col);

    % Partial pivoting
    [max_val, max_row] = max(abs(M(col:n, col)));
    max_row = max_row + col - 1;
    if max_row ~= col
      fprintf('Pivoteamento: trocando linha %d com linha %d\n', col, max_row);
      M([col max_row], :) = M([max_row col], :);
      disp('Matriz após pivoteamento:');
      disp(M);
    end

    % Normalization
    pivot = M(col, col);
    fprintf('Normalização: dividindo linha %d por pivot = %.4f\n\n', col, pivot);
    M(col, :) = M(col, :) / pivot;
    disp('Matriz após normalização:');
    disp(M);

    % Elimination
    for row = 1:n
      if row ~= col
        multiplier = M(row, col);
        fprintf('\n');
        M(row, :) = M(row, :) - multiplier * M(col, :);
        disp('Matriz após eliminação:');
        disp(M);
      end
    end
    disp('----------------------------------');
  end

  % Extract results
  A_inv = M(:, n+1:2*n);
  b_solution = M(:, 2*n+1);
end

A = [-130  30  0;
       90 -90  0;
       40  60 -120];

b = [-200;
     0;
     -500];

[A_inv, b_solution] = gauss_jordan(A, b);

disp('==================================');
disp('Resultados Finais:');

disp('Matriz inversa:');
disp(A_inv);
disp('Solução do sistema:');
fprintf('   x1 = %.0f\n', b_solution(1));
fprintf('   x2 = %.0f\n', b_solution(2));
fprintf('   x3 = %.0f\n', b_solution(3));
