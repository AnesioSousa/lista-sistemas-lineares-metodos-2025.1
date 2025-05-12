% Passo 1: Definir os valores de x
x1 = 4;
x2 = 2;
x3 = 7;

% Passo 2: Construir a matriz de Vandermonde
V = [x1^2, x1, 1;
     x2^2, x2, 1;
     x3^2, x3, 1];

disp('Matriz V (Vandermonde):');
disp(V);

% Função de Eliminação de Gauss-Jordan
function [V_inv] = gauss_jordan(A)
  [n, m] = size(A);
  if n ~= m
    error('A matriz deve ser quadrada');
  end

  M = [A, eye(n)];  % Matriz aumentada [A|I]
  disp('Matriz Inicial Aumentada [A|I]:');
  disp(M);

  % Eliminação de Gauss-Jordan
  for col = 1:n
    fprintf('\n=== Coluna %d ===\n', col);

    % Pivoteamento parcial
    [max_val, max_row] = max(abs(M(col:n, col)));
    max_row = max_row + col - 1;
    if max_row ~= col
      fprintf('Pivoteamento: trocando linha %d com linha %d\n', col, max_row);
      M([col max_row], :) = M([max_row col], :);
      disp('Matriz após pivoteamento:');
      disp(M);
    end

    % Normalização
    pivot = M(col, col);
    fprintf('Normalização: dividindo linha %d por pivô = %.4f\n\n', col, pivot);
    M(col, :) = M(col, :) / pivot;
    disp('Matriz após normalização:');
    disp(M);

    % Eliminação
    for row = 1:n
      if row ~= col
        multiplier = M(row, col);
        fprintf('\nEliminação na linha %d: subtraindo %.4f * linha %d\n', row, multiplier, col);
        M(row, :) = M(row, :) - multiplier * M(col, :);
        disp('Matriz após eliminação:');
        disp(M);
      end
    end
    disp('----------------------------------');
  end

  % Extrair a matriz inversa
  V_inv = M(:, n+1:end);  % Parte direita é a inversa de A
end

% Passo 3: Somar os valores absolutos por linha (soma linha a linha)
somas_linhas_V = sum(V_abs, 2);
disp('Soma dos valores absolutos por linha de V:');
disp(somas_linhas_V);

% Passo 4: Determinar a norma infinito de V (máximo entre as somas de linha)
norma_V = max(somas_linhas_V);
disp(['Norma infinito de V: ', num2str(norma_V)]);

% Passo 5: Calcular a inversa da matriz V usando Gauss-Jordan
V_inv = gauss_jordan(V);
disp('Matriz inversa de V (V^{-1}):');
disp(V_inv);

% Passo 6: Calcular a norma infinito da inversa de V
Vinv_abs = abs(V_inv);
disp('Matriz |V^{-1}| (valores absolutos da inversa):');
disp(Vinv_abs);

% Passo 7: Somar os valores absolutos por linha da inversa
somas_linhas_Vinv = sum(Vinv_abs, 2);
disp('Soma dos valores absolutos por linha de V^{-1}:');
disp(somas_linhas_Vinv);

% Passo 8: Determinar a norma infinito da inversa
norma_Vinv = max(somas_linhas_Vinv);
disp(['Norma infinito de V^{-1}: ', num2str(norma_Vinv)]);

% Passo 9: Calcular o número de condição pela norma infinito
cond_infinito = norma_V * norma_Vinv;
disp(['Número de condição (usando norma infinito): ', num2str(cond_infinito)]);

