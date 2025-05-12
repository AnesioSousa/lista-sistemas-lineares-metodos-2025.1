clc
clear all
close all

function x = gaussElimination(A, b)
    n = length(b);
    Ab = [A b]; % Matriz aumentada

    % Exibir matriz inicial
    disp('Matriz aumentada inicial:');
    printAugmentedMatrix(Ab);

    % Eliminação de Gauss
    for k = 1:n-1
        for i = k+1:n
            fator = Ab(i,k) / Ab(k,k);
            Ab(i,k:n+1) = Ab(i,k:n+1) - fator * Ab(k,k:n+1);
        end
        disp(['Matriz aumentada após eliminação na etapa ', num2str(k), ':']);
        printAugmentedMatrix(Ab);
    end

    % Substituição Recursiva
    x = zeros(n,1);
    x(n) = Ab(n,n+1) / Ab(n,n);
    fprintf('x%d = %.4f\n', n, x(n));
    for i = n-1:-1:1
        soma = Ab(i,n+1);
        for j = i+1:n
            soma = soma - Ab(i,j) * x(j);
        end
        x(i) = soma / Ab(i,i);
        fprintf('x%d = %.4f\n', i, x(i));
    end
end

function printAugmentedMatrix(Ab)
    [n, m] = size(Ab);
    for i = 1:n
        for j = 1:m
            fprintf('%10.4f ', Ab(i,j));  % Alinha as colunas com 4 casas decimais
        end
        fprintf('\n');
    end
    fprintf('\n');
end

% Sistema de equações
A = [0.52 0.2 0.25; 0.3 0.5 0.2; 0.18 0.3 0.55];
b = [4800; 5800; 5700];

% Resolver o sistema sem pivotamento
x = gaussElimination(A, b);

% Exibir a solução final
disp('Solução final:');
disp(x);

