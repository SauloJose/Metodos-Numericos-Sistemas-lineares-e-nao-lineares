%Método dos gradientes conjugados
function [Xr it]=GradConjSist(A,B,Xo,erro)
  % A: Matriz dos coeficientes do sistemas
  % B: Matriz dos resultados do sistema;
  % Erro: Erro desejado para a aproximação do resultados
  % Xo: Matriz solução com aproximação inicial.
  % Xr: Matriz solução aproximada
  % it: Número de iterações necessárias para a aproximação
  %Criado por: Saulo José Almeida Silva.

  %Obs: Esse método só funciona com matrizes certas.

  %Iniciando variáveis uteis
  it=1;
  %tamanho da matriz A.
  [m n]=size(A);
  if m~=n, error('A matriz deve ser quadrada');endif

  X=[Xo];
  r=[B-A*Xo]; %Gradiente
  d=[B-A*Xo]; %direção do passo negativo.
  b=[0]; %modificador do d:
  t=[0];%passo


  if (B-A*Xo)==zeros(n,1)
    error('O valor X0 já é solução!');
  endif

  %looping de iteração
  do

    fprintf("Iteração [%d]\n",it);
      b(it)=mult(r(:,it),r(:,it))/mult(A*d(:,it),d(:,it))
      X(:,it+1)=X(:,it)+b(it)*d(:,it)%próximo valor
      r(:,it+1)=r(:,it)-b(it)*A*d(:,it)
      b(it+1)=mult(r(:,it+1),r(:,it+1))/mult(r(:,it),r(:,it))
      d(:,it+1)=r(:,it+1)+b(it+1)*d(:,it)

      it=it+1;

  until it>= 200 || norm(r(:,it-1),2)<erro
  it=it-1;
  Xr=X(:,it);
  it=it-1;

endfunction

