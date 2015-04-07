function x = linprog_gen(f,A,b,B,c)
  n = size(A,2);
  m = size(A,1);
  AS = bsxfun(@times,sign(b),[A,eye(m)]);
  bS = sign(b).*b;
  n = size(A,2);
  P = sparse( ...
    [repmat(1:n,1,2) n+(1:m)], ...
    [1:2*n 2*n+(1:m)], ...
    [ones(1,n),-ones(1,n),ones(1,m)], ...
    n+m,2*n+m);
  ASP = AS*P;
  fSP = [P(1:n,1:2*n)'*f;ones(m,1)];
  cc = fSP;
  if isempty(B)
    BSP = [];
  else
    BS = [B zeros(size(B,1),n)];
    BSP = BS*P;
  end
  AA = [ASP BSP];
  bb = [bS c];
  %xxs = linprog(cc,[],[],AA,bb,zeros(2*n+m,1));
  xxs = linprog_gs(AA,bb,cc,0);
  x = P(1:n,:)*xxs;
end
