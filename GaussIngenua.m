function [x detA] =GaussIngenua(A,B)
  % A: Matriz dos coeficientes do sistemas
  % B: Matriz dos resultados do sistema;
  % X: Matriz solução aproximada
  %Criado por: Saulo José Almeida Silva.

  %Tamanho da matriz
  [m,n]=size(A);
  if m~=n, error('A matriz deve ser quadrada');end
  np=n+1;
  Aum=[A B];
  pivot=0;
  detA=1;

  %Eliminação progressiva
  for k=1:n-1
    %Revezando
    for i=k+1:n
      fator=Aum(i,k)/Aum(k,k);
      Aum(i,k:np)=Aum(i,k:np)-fator*Aum(k,k:np);
    endfor
  endfor
  Aum

  %substituição regressiva.
  x=zeros(n,1);
  x(n)=Aum(n,np)/Aum(n,n); %bn/ancestor

  for i=n-1:-1:1
    x(i)=(Aum(i,np)-Aum(i,i+1:n)*x(i+1:n))/Aum(i,i);
  endfor

  %calculando o determinante por gauss
  for i=1:n
    detA=detA*Aum(i,i)*(-1)^pivot;
  endfor

endfunction
