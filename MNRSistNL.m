%Método de Newton-RAFAEL
function [X it Fr dFr X1]=MNRSistNL(X0,erro)
  %[Xr it] = AproxSucessivas(A,B,Xo,erro)
  % Entradas:
  %Xo: Vetor de aprximação inicila
  %erro: erro esperado para a resposta
  %Saidas:
  %Xr: Vetor solução
  %it: Número de iterações necessárias
  %Por: Saulo José Almeida Silva

  %vetor solução
  Xr=[X0];
  Dxr=[X0];

  %definindo funções de iteração
  nx=@(x,y,z) x^2+y^2+z^2-9;
  ny=@(x,y,z) x*y*z-1;
  nz=@(x,y,z) x+y-z^2;

  %definindo derivadas parciais
  %função em x
  dnxx=@(x,y,z) 2*x;
  dnxy=@(x,y,z) 2*y;
  dnxz=@(x,y,z) 2*z;

  %função em y
  dnyx=@(x,y,z) y*z;
  dnyy=@(x,y,z) x*z;
  dnyz=@(x,y,z) x*y;


  %função em z
  dnzx=@(x,y,z) 1;
  dnzy=@(x,y,z) 1;
  dnzz=@(x,y,z) -2*z;

  %definindo função F
  F=@(x,y,z) [nx(x,y,z);ny(x,y,z);nz(x,y,z)];

  %Definindo função derivada.
  dF=@(x,y,z) [dnxx(x,y,z) dnxy(x,y,z) dnxz(x,y,z);dnyx(x,y,z) dnyy(x,y,z) dnyz(x,y,z);dnzx(x,y,z) dnzy(x,y,z) dnzz(x,y,z)];

  %inversa da função derivada
  invdF=@(x,y,z) inv(dF(x,y,z));

  %variavel de iteração
  it=1;
  U=zeros(3,3);
  do
     it=it+1;
     Xr(:,it)=Xr(:,it-1)-invdF(Xr(1,it-1),Xr(2,it-1),Xr(3,it-1))*F(Xr(1,it-1),Xr(2,it-1),Xr(3,it-1));
     Dxr(:,it)=Xr(:,it)-Xr(:,it-1);
  until norm(Dxr(:,it),2)<erro || it >=30

  X=Xr(:,it);
  %Imprimindo primeira iteração.
  Fr=F(Xr(1,1),Xr(2,1),Xr(3,1));
  dFr=dF(Xr(1,1),Xr(2,1),Xr(3,1));
  X1=Xr(:,2);

  %valores finais da função
  fprintf("Valores da função\n");
  nx(Xr(1,it),Xr(2,it),Xr(3,it))
  ny(Xr(1,it),Xr(2,it),Xr(3,it))
  nz(Xr(1,it),Xr(2,it),Xr(3,it))
  fprintf("\nResultados:");

endfunction

