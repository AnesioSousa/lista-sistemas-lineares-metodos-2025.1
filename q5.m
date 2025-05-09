clc
clear all
close all

% Função de eliminação de Gauss com prints e registro de erros
function [X, errors] = gauss_elimination_with_convergence(A, B)
    [n, m] = size(A);

    % Montando a matriz aumentada
    AugmentedMatrix = [A, B];
    disp('Matriz aumentada inicial:');
    disp(AugmentedMatrix);

    % Inicializando vetor de erros
    errors = [];

    % Fase de eliminação
    for i = 1:n
        for j = i+1:n
            fator = AugmentedMatrix(j, i) / AugmentedMatrix(i, i);
            AugmentedMatrix(j, :) = AugmentedMatrix(j, :) - fator * AugmentedMatrix(i, :);
            fprintf('Matriz aumentada após zerar o elemento (%d, %d):\n', j, i);
            disp(AugmentedMatrix);
        end
    end

    % Fase de retro-substituição
    X = zeros(n, 1);
    for i = n:-1:1
        X(i) = (AugmentedMatrix(i, end) - AugmentedMatrix(i, 1:end-1) * X) / AugmentedMatrix(i, i);
        fprintf('Solução parcial após calcular x_%d:\n', i);
        disp(X);

        % Calcular erro relativo e armazenar
        error = norm(A * X - B) / norm(B);
        errors = [errors; error];
    end
end

% Definindo a matriz dos coeficientes A e o vetor B
A = [10 2 -1;
     -3 -6 2;
     1 1 5];

B = [27;
     -61.5;
     -21.5];

% Resolvendo o sistema de equações e registrando os erros
[X, errors] = gauss_elimination_with_convergence(A, B);

% Exibindo o resultado final
disp('A solução do sistema de equações é:')
disp(X)

% Verificando o resultado
A_inv = inv(A);
I = A * A_inv;
disp('A matriz identidade resultante da multiplicação de A por sua inversa é:')
disp(I)

% Plotando o gráfico de convergência
figure;
plot(errors, '-o');
xlabel('Iteração');
ylabel('Erro Relativo');
title('Gráfico de Convergência da Eliminação de Gauss');
grid on;