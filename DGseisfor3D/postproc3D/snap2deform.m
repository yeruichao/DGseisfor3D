function snap2deform(snapname,outname,eid,src,srcr,amp)

dir='/scratch/conte/r/rye/snap/'
outdir='/scratch/conte/r/rye/snap/'
fname=[dir,snapname,'_Para.txt'];
Para=load(fname);
Nbele =Para(1);
Np    =Para(2);
Nsele =Para(3);
Nproc =Para(4);
Nsnap =Para(5);
Nnod  = Nbele*Np;
Nele  = Nbele*Nsele;
DNele = zeros(Nproc,1);

for i=1:Nproc
    fidb=fopen([dir,snapname,'_globID',num2str(i-1),'.txt'],'r');
    DNele(i)=fscanf(fidb,'%d',1);
    fclose(fidb);
end
if(sum(DNele)~=Nbele)
    disp('Something wrong with the number of elements.');
end

p=zeros(Nnod,3);offset=0;
for i=1:Nproc
    fidb=fopen([dir,snapname,'_Coord',num2str(i-1),'.dat'],'r');
    p(offset+1:offset+DNele(i)*Np,1)=fread(fidb,DNele(i)*Np,'float');
    p(offset+1:offset+DNele(i)*Np,2)=fread(fidb,DNele(i)*Np,'float');
    p(offset+1:offset+DNele(i)*Np,3)=fread(fidb,DNele(i)*Np,'float');
    offset=offset+DNele(i)*Np;
    fclose(fidb);
end

disp('Coordinate done')

t=zeros(4*Nele,1);offset=0;counter=0;
for i=1:Nproc
    fidb=fopen([dir,snapname,'_Conn',num2str(i-1),'.txt'],'r');
    t(offset+1:offset+DNele(i)*Nsele*4)=fscanf(fidb,'%d',DNele(i)*Nsele*4)+counter;
    offset=offset+DNele(i)*Nsele*4;counter=counter+DNele(i)*Np;
    fclose(fidb);
end
t=reshape(t,4,Nele);

disp('Connection done')
t=t';

id=false(Nsele,Nbele);
id(:,eid)=true;
id=id(:);
[t,nid]=pickless(p,t,id);
p=p(nid,:);
nnod=size(p,1);
nele=size(t,1);

r=sqrt((p(:,1)-src(1)).^2+(p(:,2)-src(2)).^2+(p(:,3)-src(3)).^2);
w=ones(length(r),1);
w(r<=srcr)=0;
id=r>srcr & r<srcr*2;
w(id)=(sin((r(id)/srcr-1.5)*pi)+1)/2;

Nnod
nnod
Nele
nele

for i=1:Nsnap
    v=zeros(Nnod,3);offset=0;
    for j=1:Nproc
      fidb=fopen([dir,snapname,'_V1_',num2str(i),'_',num2str(j-1),'.dat'],'r');
        v(offset+1:offset+DNele(j)*Np,1)=fread(fidb,DNele(j)*Np,'float');
        fclose(fidb);
      fidb=fopen([dir,snapname,'_V2_',num2str(i),'_',num2str(j-1),'.dat'],'r');
        v(offset+1:offset+DNele(j)*Np,2)=fread(fidb,DNele(j)*Np,'float');
        fclose(fidb);
      fidb=fopen([dir,snapname,'_V3_',num2str(i),'_',num2str(j-1),'.dat'],'r');
        v(offset+1:offset+DNele(j)*Np,3)=fread(fidb,DNele(j)*Np,'float');
        fclose(fidb);
        offset=offset+DNele(j)*Np;
    end
    max(v)
    v=v(nid,:);
    v=v.*(w*ones(1,3));
    v=p+amp*v;
    vtk_write([dir,snapname,num2str(i),'.vtk'],'snap',v,t);
    clear v;
    disp([num2str(i),'th snapshot of vector done.'])
end

vtuxml_end(fid);

end
