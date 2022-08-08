  #função para método LU
  function [L U x]=DecLU(A,B)
    %   Função para o método LU
    % A : Matriz dos Coeficientes
    % B : Matriz dos Resultados
    % L : Matriz triangular INFERIOR
    % U : Matriz triangular SUPERIOR
    % X : Matriz solução do sistema.
    % Criado por: Saulo José

    %Tamanho da matriz
    [m n]=size(A);
    if m~=n, error('A matriz deve ser quadrada');end

    %gerando matrizes U e L para manipular.
    U=A;
    Aum=[A B];
    L=eye(n);
    d=zeros(n,1);
    x=zeros(n,1);

    %PARTE 1: manipulação de A para chegar em U e L (OK)
    for k=1:n-1
        %pivoteamento => evita divisões por zeros
        [maior,i]=max(abs(U(k:n,k)));
        ipr=i+k-1;
        if ipr~=k
        U([k,ipr],:)=U([ipr,k],:);
        endif

     %Gerando matrizes U e V
      for i=k+1:n
        fator = U(i,k)/U(k,k); %GERANDO fatores
        L(i,k)=fator;% Gerando L
        U(i,k:n)=U(i,k:n)-fator*U(k,k:n);%Gerando u
      endfor
    endfor

    %PARTE 2: Encontrar a matriz intermediária d => Ld=b =>
    B = (L*U/A)*B

     d(1)=B(1);
     for i=2:n
        d(i)=B(i)-L(i,1:i-1)*d(1:i-1)
     endfor

    %Parte 3: Resolver U*x=d
     x(n)=d(n)/U(n,n)

     for i=n-1:-1:1
       x(i)=(d(i)-U(i,i+1:n)*x(i+1:n))/U(i,i)
     endfor
endfunction

