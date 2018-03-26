function [t1,nid]=pickless(p,t,eid)

if(size(p,1)==3)
  p=p';
end
if(size(t,1)==4)
  t=t';
end

t1=t(eid,:);
Nele=size(t1,1);

Nnode=size(p,1);
flag=zeros(Nnode,1);
t1=t1(:);
flag(t1)=1;
nid=find(flag==1);
Nnode=size(nid,1);
flag(nid(:))=1:Nnode;
t1=flag(t1);
t1=reshape(t1,Nele,4);

end
