clc
clear all
close all

function [x, det] = gaussEliminationWithPivot(A, b)
    sinal = 1; % Para acompanhar as trocas de linha
    n = length(b);
    Ab = [A b]; % Matriz aumentada

    % Exibir matriz inicial
    disp('Matriz aumentada inicial:');
    disp(Ab);

    % Eliminação com pivotamento parcial
    for k = 1:n-1
        % Pivoteamento Parcial
        [maior, I] = max(abs(Ab(k:n, k)));
        ipr = I + k - 1;
        if ipr ~= k
            sinal = -sinal; % Cada troca de linha inverte o sinal do determinante
            % Troca de linhas
            Aux = Ab(k,:);
            Ab(k,:) = Ab(ipr,:);
            Ab(ipr,:) = Aux;
            disp(['Troca de linha ', num2str(k), ' com linha ', num2str(ipr)]);
        end

        % Exibir matriz após pivoteamento
        disp(['Matriz aumentada após pivoteamento da coluna ', num2str(k), ':']);
        disp(Ab);

        % Eliminação
        for i = k+1:n
            fator = Ab(i,k) / Ab(k,k);
            Ab(i,k:n+1) = Ab(i,k:n+1) - fator * Ab(k,k:n+1);
        end

        % Exibir matriz após eliminação
        disp(['Matriz aumentada após eliminação na etapa ', num2str(k), ':']);
        disp(Ab);
    end

    det = sinal * prod(diag(Ab(:, 1:n)));

    % Substituição Recursiva
    x = zeros(n,1);
    x(n) = Ab(n,n+1) / Ab(n,n);
    disp(['x', num2str(n), ' = ', num2str(x(n))]); % Exibir o valor da última variável
    for i = n-1:-1:1
        soma = Ab(i,n+1);
        for j = i+1:n
            soma = soma - Ab(i,j) * x(j);
        end
        x(i) = soma / Ab(i,i);
        disp(['x', num2str(i), ' = ', num2str(x(i))]); % Exibir o valor de cada variável
    end
end

% Sistema de equações
A = [2 -6 -1; -3 -1 7; -8 1 -2];
b = [-38; -34; -20];
% Resolver o sistema com pivotamento
[x, det] = gaussEliminationWithPivot(A, b);

% Exibir a solução final
disp('Solução final:');
disp(x);

disp(['Determinante: ', num2str(det)]);
