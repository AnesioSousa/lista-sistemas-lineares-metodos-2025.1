clc;
clear all;
close all;

% Matriz original
A = [8  2 -10;
    -9  1   3;
    15 -1   6];

% Passo 1: Normalizar a matriz A
n = size(A, 1); % Número de linhas
Anorm = A; % Inicializando a matriz normalizada

% Construindo a normalização de cada linha como um sistema linear
% Método: Normalização direta das linhas
disp('--- Início da Normalização da Matriz A ---');
maxElems = zeros(n, 1); % Vetor dos maiores elementos das linhas
for i = 1:n
    maxElems(i) = max(abs(A(i, :))); % Maior elemento da linha i
    Anorm(i, :) = A(i, :) / maxElems(i);
end

% Resultado da matriz normalizada
disp('Matriz normalizada:');
disp(Anorm);

% Passo 2: Calcular as normas da matriz normalizada
norma1 = max(sum(abs(Anorm), 1));
normaInfinito = max(sum(abs(Anorm), 2));

% Exibindo as normas calculadas
fprintf('Norma 1 (|A|_1): %.4f\n', norma1);
fprintf('Norma infinito (|A|_inf): %.4f\n', normaInfinito);

% Passo 3: Implementar o método de Gauss-Seidel para convergência
disp('--- Início do Método de Gauss-Seidel ---');
tol = 1e-6; % Tolerância de erro
maxIter = 100; % Máximo de iterações
x = zeros(n, 1); % Vetor inicial de solução
b = maxElems; % Vetor do maior elemento de cada linha (normalização)
B = diag(maxElems); % Matriz diagonal (escalares)
C = Anorm; % Matriz após normalização
iterations = []; % Vetor para armazenar o erro em cada iteração

% Método de Gauss-Seidel
for k = 1:maxIter
    x_prev = x; % Solução da iteração anterior
    for i = 1:n
        sum1 = 0;
        for j = 1:n
            if j ~= i
                sum1 = sum1 + C(i, j) * x(j);
            end
        end
        x(i) = (b(i) - sum1) / B(i, i);
    end
    error = norm(x - x_prev, Inf) / norm(x, Inf);
    iterations = [iterations, error];
    if error < tol
        break; % Convergência atingida
    end
end

% Exibindo a solução final e número de iterações
disp('Solução final (fatores de normalização):');
disp(x);
fprintf('Número total de iterações: %d\n', k);

% Passo 4: Traçar o gráfico de convergência
figure;
plot(1:length(iterations), iterations, '-o', 'LineWidth', 2);
xlabel('Iteração');
ylabel('Erro Relativo');
title('Convergência do Método de Gauss-Seidel');
grid on;