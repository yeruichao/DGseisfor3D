function showsnap3D(snapname,meshname,outname,varname,i_snap)

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

Scalname=cell(Nsnap,1);
Vectname=cell(Nsnap,1);

Nsnap=5
N=2;

for i=1:Nproc
    fidb=fopen([dir,snapname,'_globID',num2str(i-1),'.txt'],'r');
    DNele(i)=fscanf(fidb,'%d',1);
    fclose(fidb);
end
if(sum(DNele)~=Nbele)
    disp('Something wrong with the number of elements.');
end

i=i_snap;
    v=zeros(Nnod,1);offset=0;
    for j=1:Nproc
      fidb=fopen([dir,snapname,'_',varname,'_',num2str(i),'_',num2str(j-1),'.dat'],'r');
        v(offset+1:offset+DNele(j)*Np)=fread(fidb,DNele(j)*Np,'float');
        offset=offset+DNele(j)*Np;
        fclose(fidb);
    end
    max(v)
    disp([num2str(i),'th snapshot of scaler done.'])
maxv=max(v)
minv=min(v)
v=reshape(v,Np,Nbele)';
v=v(:,[1,N+1,(N+1)*(N+2)/2,Np]);

[p,t]=read_mesh3d(meshname);
v1=zeros(size(p,1),1);
v1(t)=v;

maxv=max(v1)
minv=min(v1)

%[t1,eid,nid,nnode,nnele]=pickeles(p,t);
eid=v(:,1)>1 | v(:,1)<-1;
[t1,nid]=pickless(p,t,eid);
fname=[dir,snapname,'_',varname,'_',num2str(i_snap),'.vtk']
vtk_write(fname,varname,p(nid,:),t1,v(eid,:),1);

end
