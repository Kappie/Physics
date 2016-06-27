function [u,s,v]=tensorsvd(T,ll,rl,n,absorb,cut)
%[u,s,v]=tensorsvd(T,ll,rl,n,absorb,cut)
%do svd for a tensor
%T: tensor
%ll: left legs
%rl: right legs
%n: numer of singular values: use inf for all
%absorb: 'l': left absorb
%'m': split absorb (sqrt of singular values)
%'r': right absorb
%'b': absorb singular values on both sides
%'n': no absorption of singular values
%cut: cutoff for smallest SVDs
%structure of output:
%order of legs u: leftlegs, single leg
%order of legs v: rightlegs, single leg
%CONVENTION: the conjugate of v is taken such that the original tensor is obtained
%by contracting u.s.v without taking transpose/conjugate

[Tm,dl,dr]=lreshape(T,ll,rl);

%default: no absorption
if(nargin<5)
    absorb='n';
end

%default cutoff
if(nargin<6)
    cut=1e-14;
end

%Default:
[u,s,v]=svdsabsorb(Tm,n,absorb,cut);

ss=size(u);
u=reshape(u,[dl,ss(2)]);
v=reshape(v,[dr,ss(2)]); 

v=conj(v);

if(nargout==2)
    s=v;
end

end