#Método dos Gradientes puros para sistemas lineares
function [Y it]=GradPSist(A,B,Xi,erro)
  % A: Matriz dos coeficientes do sistemas
  % B: Matriz dos resultados do sistema;
  % Erro: Erro desejado para a aproximação do resultados
  % Xi: Matriz solução com aproximação inicial.
  % X: Matriz solução aproximada
  % it: Número de iterações necessárias para a aproximação
  %Criado por: Saulo José Almeida Silva.

  %Iniciando variáveis uteis
  it=1;
  %tamanho da matriz A.
  [m n]=size(A);
  if m~=n, error('A matriz deve ser quadrada');end

  X=[Xi];
  g=[zeros(n,1)]; %Gradiente
  d=[zeros(n,1)]; %direção do passo negativo.
  t=0;%passo

  g(:,1)= A*X(:,1)-B;
  if g(:,1)==zeros(n,1)
    error('O valor X0 já é solução!');
  endif

  %looping de iteração
  do
    it=it+1;
    g(:,it-1)=A*X(:,it-1)-B;
    d(:,it-1)=-g(:,it-1);
    t(it-1)=mult(g(:,it-1),g(:,it-1))/mult(A*g(:,it-1),g(:,it-1));
    X(:,it)=X(:,it-1)+t(it-1)*d(:,it-1);

until it>= 100 || norm(g(:,it-1),2)<erro
  %Primeira iteração
  fprintf("Primeira iteração\n")
  fprintf("X(0)");
  X(:,1)
  fprintf("Gradiente(1)");
  g(:,1)
  fprintf("-Grad(1)");
  d(:,1)
  fprintf("passo(1)");
  t(1)
  fprintf("Valor da primeira iteração");
  X(:,2)
  Y=X(:,it);
endfunction
