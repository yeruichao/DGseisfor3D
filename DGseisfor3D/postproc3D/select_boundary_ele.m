function Tid=select_boundary_ele(N)

p=[0,0,0;1,0,0;0,1,0;0,0,1];
tet=[1,2,3,4];
nodal=load('porder2nodal.mat','porder');
r=nodal.porder{N}.r;
s=nodal.porder{N}.s;
t=nodal.porder{N}.t;
x = 0.5*(-(1+r+s+t)*p(tet(:,1),1)'...
             +(1+r)*p(tet(:,2),1)'...
             +(1+s)*p(tet(:,3),1)'...
             +(1+t)*p(tet(:,4),1)');
y = 0.5*(-(1+r+s+t)*p(tet(:,1),2)'...
             +(1+r)*p(tet(:,2),2)'...
             +(1+s)*p(tet(:,3),2)'...
             +(1+t)*p(tet(:,4),2)');
z = 0.5*(-(1+r+s+t)*p(tet(:,1),3)'...
             +(1+r)*p(tet(:,2),3)'...
             +(1+s)*p(tet(:,3),3)'...
             +(1+t)*p(tet(:,4),3)');

tt=tet_subelement(N);
Tid=zeros(size(tt,1),4);

id=x<1e-3;
tid=zeros(size(tt));
tid(id(tt))=1;
tid=sum(tid,2);
tid=tid==3;
Tid(tid,2)=1;

id=y<1e-3;
tid=zeros(size(tt));
tid(id(tt))=1;
tid=sum(tid,2);
tid=tid==3;
Tid(tid,3)=1;

id=z<1e-3;
tid=zeros(size(tt));
tid(id(tt))=1;
tid=sum(tid,2);
tid=tid==3;
Tid(tid,4)=1;

id=x+y+z>1-1e-3;
tid=zeros(size(tt));
tid(id(tt))=1;
tid=sum(tid,2);
tid=tid==3;
Tid(tid,1)=1;

