function [u,s,v]=svdsabsorb(T,n,mode,cut)
%T input matrix
%n: number of singular values (set inf for all)
%mode: 'l': left absorb
%'m': split absorb
%'r': right absorb
%'b': absorb on both sides
%'n': no absorption of singular values

%set precision for cutoff used in svd
%global SVDPRECCUT 

%default cutoff
if(nargin<4)
    cut=1e-14;
end
 
%if(~isempty(SVDPRECCUT))
%    cut=SVDPRECCUT;
%end

% DO SVD   
try
    [u,s,v]=svd(T,'econ');
catch me
    disp(me.message);
    disp('ERROR OCCURED IN SVD, TRY SVDS instead:')
    [u,s,v]=svds(T,n);
end
       
n=min(n,size(s,1));

%check cutoff: truncate very small singular values
ds=diag(s);
ds=ds(1:n);
icut=find(ds/max(ds)<cut,1);
if(~isempty(icut))
    n=icut-1;
end

u=u(:,1:n);
v=v(:,1:n);
s=s(1:n,1:n);

% absorb singular values?
if(isequal(mode,'m'))
   sq=diag(sqrt(diag(s)));
   u=u*sq;
   v=v*sq;
elseif(isequal(mode,'l'))
   u=u*s;
elseif(isequal(mode,'b'))
   u=u*s;
   v=v*s;    
elseif(isequal(mode,'r'))    
   v=v*s;
elseif(isequal(mode,'n'))    
   %nothing 
else
   error('wrong mode');
end

end


