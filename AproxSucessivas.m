%Aproximações sucessivas e newton rapson
function [Xr it] = AproxSucessivas(Xo,erro)
  %[Xr it] = AproxSucessivas(A,B,Xo,erro)
  % Entradas:
  %Xo: Vetor de aprximação inicila
  %erro: erro esperado para a resposta
  %Saidas:
  %Xr: Vetor solução
  %it: Número de iterações necessárias
  %Por: Saulo José Almeida Silva

  %inicializando variavel de iteração
  it=1;

  %funções de iteração (Aqui você coloca suas funções)
  f1=@(x,y) (-2-2*y)^(1/2);
  f2=@(x,y) -((21-x^3)/3)^(1/2);
  Fop=@(x,y) [f1(x,y);f2(x,y)];

  %Matriz para salvar os valores das iterações
  Xr=[Xo];
  do
    it=it+1;
    Xr(:,it)=Fop(Xr(1,it-1),Xr(2,it-1))
  until abs(Xr(1,it)-Xr(1,it-1))<erro && abs(Xr(2,it)-Xr(2,it-1))<erro || it>20

endfunction

