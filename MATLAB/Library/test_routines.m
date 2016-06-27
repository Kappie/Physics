a=rand(2,3,5,2,1);

%reshape into matrix
am=lreshape(a,[1 3],[2 4 5]);

%do a tensorsvd
[u,s,v] = tensorsvd(a,[1 3],[2 4 5],inf,'n');

%check result:
res = ncon({u,s,v},{[-1 -3 1],[1 2],[-2 -4 -5 2]});
check=a-res;
error1=sum(abs(check(1:end)))

%do a tensorsvd with truncation:
[u,s,v] = tensorsvd(a,[1 3],[2 4 5],5,'n');

%check result:
res = ncon({u,s,v},{[-1 -3 1],[1 2],[-2 -4 -5 2]});
check=a-res;
error2=sum(abs(check(1:end)))


%do a tensorsvd with absorbing the singular values on the left:
[u,s,v] = tensorsvd(a,[1 3],[2 4 5],6,'l');
res = ncon({u,v},{[-1 -3 1],[-2 -4 -5 1]});
check=a-res;
error3=sum(abs(check(1:end)))


%do a tensorsvd with cutoff at 0.1
[u,s,v] = tensorsvd(a,[1 3],[2 4 5],inf,'n',0.1);

%check result:
res = ncon({u,s,v},{[-1 -3 1],[1 2],[-2 -4 -5 2]});
check=a-res;
error4=sum(abs(check(1:end)))