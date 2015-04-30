function [sol,exato] = elementos_finitos_esparso(N)
% ELEMENTOS_FINITOS_ESPARSO
% 
% sol = elementos_finitos_esparso(N)
%
% Usando elementos finitos lineares por partes,
% resolve o problema de valores de contorno -u'' = f, no intervalo [0,1],
% com condicoes de contorno u(0) = u(1) = 0, para f(x) = - x^3 + 5.
% Define matriz do sistema como esparsa.
% 
% Entrada:
%   N: N+1 eh o numero de pontos para discretizar o intervalo [0,1]
% Saida:
%   sol: vetor com valores da solucao aproximada nos pontos de malha 
%       {h,2h,...,(N-1)h}
%   exato: vetor com solucao exata 

% define espacamento da malha
h = 1/N;

% define lado direito
i = (1:(N-1))';
b = (h^4)*((i.^5)/10-((i-1).^5)/20-((i+1).^5)/20)+...
    h*(5*(i.^2)-(5/2)*(i-1).^2-(5/2)*(i+1).^2+10);

% define matriz do sistema linear
I = [2:(N-1) 1:(N-1) 1:(N-2)];
J = [1:(N-2) 1:(N-1) 2:(N-1)];
Val = [-(1/h)*ones(1,N-2) (2/h)*ones(1,N-1) -(1/h)*ones(1,N-2)];
A = sparse(I,J,Val,N-1,N-1,3*(N-1));

% resolve sistema linear
sol = A\b;

% retorna tambem solucao exata do problema
x = i*h;
exato = (x.^5)/20-(5/2)*(x.^2) + (49/20)*x;