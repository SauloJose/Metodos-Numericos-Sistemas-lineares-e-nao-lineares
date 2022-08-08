#Decomposição Tridiag
function x=Tridiag(e,f,g,r)
  %Tridiag: Solução de sistemas de equação tridiagonais
  % x=Tridiag(e,f,g,r): soluções de sistemas de equações tridimensionais
  %entreda:
  % e=vetor subdiagonal
  % f=vetor diagonal
  % g= vetor superdiagonal
  % r=vetor do lado direito
  %saida:
  % x = vetor solução

  n=length(f);
  %eliminação progressiva;
  for k=2:n
    factor=e(k)*f(k-1);
    f(k)=f(k)-factor*g(k-1);
    r(k)=r(k)-factor*r(k-1);
  endfor

  x(n)=r(n)/f(n);
  for k=n-1:-1:1
    x(k)=(r(k)-g(k)*x(k+1))/f(k);
  endfor

endfunction
