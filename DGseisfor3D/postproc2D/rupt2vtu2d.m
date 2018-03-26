function rupt2vtu2d(snapname,outfile,Np,Nproc,startsnap,finalsnap,intv,amp)

tt=trisubele(Np);
Nsele=size(tt,1);
Nfp=(Np+1)*(Np+2)/2;
pNp=(Np+1)*(Np+2)*(Np*3)/6;
tt=tt';

DNele = zeros(Nproc,1);
DhNele = zeros(Nproc,1);

for i=1:Nproc
    fname=[snapname,'_Rupture_',num2str(i-1,'%d'),'.dat'];
    if(not(exist(fname,'file')))
        disp(fname)
        return
    end
    fidb=fopen(fname,'r');
    globNface=fread(fidb,1,'int');
    DNele(i)=fread(fidb,1,'int');
    DhNele(i)=fread(fidb,1,'int');
    fclose(fidb);
end
%DNele
%Nfp
globNface=sum(DhNele);

Nnod  = globNface*Nfp;
Nele  = globNface*Nsele;

Nscalar_ivar=3;
Scalname_ivar=cell(Nscalar_ivar,1);
Scalname_ivar{1}='pid';
Scalname_ivar{2}='nSn_ini';
Scalname_ivar{3}='a';

Nvector_ivar=1;
Vectname_ivar=cell(Nvector_ivar,1);
Vectname_ivar{1}='tau_ini';

Nscalar_tvar=4;
Scalname_tvar=cell(Nscalar_tvar,1);
Scalname_tvar{1}='nSn';
Scalname_tvar{2}='Vt';
Scalname_tvar{3}='psi';
Scalname_tvar{4}='f';
%Scalname_tvar{5}='ruptclosed';

Nvector_tvar=3;
Vectname_tvar=cell(Nvector_tvar,1);
Vectname_tvar{1}='tauf';
Vectname_tvar{2}='tau';
Vectname_tvar{3}='dVt';

Nt=length(startsnap:intv:finalsnap);

Nscal=Nscalar_ivar+Nscalar_tvar*Nt;
Scalname=cell(Nscal,1);
for iscal=1:Nscalar_ivar
    Scalname{iscal}=Scalname_ivar{iscal};
end
for i=1:Nscalar_tvar
    for j=startsnap:intv:finalsnap
        iscal=iscal+1;
        Scalname{iscal}=[Scalname_tvar{i},'_t',num2str(j,'%04d')];
    end
end

Nvect=Nvector_ivar+Nvector_tvar*Nt;
Vectname=cell(Nvect,1);
for ivect=1:Nvector_ivar
    Vectname{ivect}=Vectname_ivar{ivect};
end
for i=1:Nvector_tvar
    for j=startsnap:intv:finalsnap
        ivect=ivect+1;
        Vectname{ivect}=[Vectname_tvar{i},'_t',num2str(j,'%04d')];
    end
end

fid=vtuxml_head_tri(outfile,Nnod,Nele,Scalname,Vectname);

p  =zeros(Nnod,3);
Sti=zeros(Nnod,3);
Sni=zeros(Nnod,1);
a  =zeros(Nnod,1);
pid=zeros(Nnod,1);
offset=0;
for i=1:Nproc
    if(DNele(i)>=0)
        pid(offset+1:offset+DhNele(i)*Nfp)=i-1.0;
        fname=[snapname,'_Rupture_',num2str(i-1,'%d'),'.dat'];
        if(not(exist(fname,'file')))
            disp(fname)
            return
        end
        fidb=fopen(fname,'r');
        fread(fidb,3,'int');
        fread(fidb,DNele(i),'int');
%        fread(fidb,DNele(i)*Nfp*2,'int');
        tmp=fread(fidb,DNele(i)*Nfp,'double');
        p(  offset+1:offset+DhNele(i)*Nfp,1)=tmp(1:DhNele(i)*Nfp);
        tmp=fread(fidb,DNele(i)*Nfp,'double');
        p(  offset+1:offset+DhNele(i)*Nfp,2)=tmp(1:DhNele(i)*Nfp);
        tmp=fread(fidb,DNele(i)*Nfp,'double');
        p(  offset+1:offset+DhNele(i)*Nfp,3)=tmp(1:DhNele(i)*Nfp);
%        tmp=fread(fidb,DNele(i)*Nfp,'double');
%        Sti(offset+1:offset+DhNele(i)*Nfp,1)=tmp(1:DhNele(i)*Nfp);
%        tmp=fread(fidb,DNele(i)*Nfp,'double');
%        Sti(offset+1:offset+DhNele(i)*Nfp,2)=tmp(1:DhNele(i)*Nfp);
%        tmp=fread(fidb,DNele(i)*Nfp,'double');
%        Sti(offset+1:offset+DhNele(i)*Nfp,3)=tmp(1:DhNele(i)*Nfp);
%        tmp=fread(fidb,DNele(i)*Nfp,'double');
%        a(  offset+1:offset+DhNele(i)*Nfp  )=tmp(1:DhNele(i)*Nfp);% vt1
%        tmp=fread(fidb,DNele(i)*Nfp,'double');
%        a(  offset+1:offset+DhNele(i)*Nfp  )=tmp(1:DhNele(i)*Nfp);% vt2
%        tmp=fread(fidb,DNele(i)*Nfp,'double');
%        a(  offset+1:offset+DhNele(i)*Nfp  )=tmp(1:DhNele(i)*Nfp);% vt3
        tmp=fread(fidb,DNele(i)*Nfp,'double');
        a(  offset+1:offset+DhNele(i)*Nfp  )=tmp(1:DhNele(i)*Nfp);% a
        tmp=fread(fidb,DNele(i)*Nfp,'double');
        Sni(offset+1:offset+DhNele(i)*Nfp  )=tmp(1:DhNele(i)*Nfp);
        fclose(fidb);
        offset=offset+DhNele(i)*Nfp;
    end
end

disp('')
disp('time invariant done')

t=Nfp*ones(3*Nsele,1)*(0:globNface-1);
t=t+tt(:)*ones(1,globNface);
t=reshape(t,3,Nsele*globNface)-1;

vtuxml_xyz(fid,p);
disp('coordinate done')
vtuxml_conn_tri(fid,t);
disp('connection done')
clear p t

St=zeros(Nnod,3,Nt);
Sf=zeros(Nnod,3,Nt);
dVt=zeros(Nnod,3,Nt);
Sn=zeros(Nnod,1,Nt);
Vt=zeros(Nnod,1,Nt);
Ut=zeros(Nnod,1,Nt);
f =zeros(Nnod,1,Nt);
%bo=zeros(Nnod,1,Nt);
offset=0;offsetV=0;
for i=1:Nproc
    if(DNele(i)>0)
        j=0;
        for s=startsnap:intv:finalsnap
            j=j+1;
            fname=[snapname,'_Rupt_',num2str(s,'%d'),'_',num2str(i-1,'%d'),'.dat'];
            if(not(exist(fname,'file')))
                disp(fname)
                return
            end
            fidb=fopen(fname,'r');
%% nSn
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            Sn( offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
%% tauf
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            Sf( offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            Sf( offset+1:offset+DhNele(i)*Nfp,2,j)=tmp(1:DhNele(i)*Nfp);
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            Sf( offset+1:offset+DhNele(i)*Nfp,3,j)=tmp(1:DhNele(i)*Nfp);
%% dVt
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            dVt(offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            dVt(offset+1:offset+DhNele(i)*Nfp,2,j)=tmp(1:DhNele(i)*Nfp);
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            dVt(offset+1:offset+DhNele(i)*Nfp,3,j)=tmp(1:DhNele(i)*Nfp);
%	    for ii=1:3
%                tmp=fread(fidb,DNele(i)*Nfp,'double');
%	    end
%            tmp=fread(fidb,DNele(i)*Nfp,'double');
%            St( offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
%            tmp=fread(fidb,DNele(i)*Nfp,'double');
%            St( offset+1:offset+DhNele(i)*Nfp,2,j)=tmp(1:DhNele(i)*Nfp);
%            tmp=fread(fidb,DNele(i)*Nfp,'double');
%            St( offset+1:offset+DhNele(i)*Nfp,3,j)=tmp(1:DhNele(i)*Nfp);
            fread(fidb,DNele(i)*Nfp,'double');
%            Vt( offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
%% f
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            f(  offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
%% psi
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            Ut( offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
%% Crack_t
            tmp=fread(fidb,DNele(i)*Nfp,'double');
            Vt( offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
%            for ii=1:24
%                fread(fidb,DNele(i)*pNp,'double');
%            end
%%% tau
%            tmp=fread(fidb,DNele(i)*Nfp,'double');
%            St( offset+1:offset+DhNele(i)*Nfp,1,j)=tmp(1:DhNele(i)*Nfp);
%            tmp=fread(fidb,DNele(i)*Nfp,'double');
%            St( offset+1:offset+DhNele(i)*Nfp,2,j)=tmp(1:DhNele(i)*Nfp);
%            tmp=fread(fidb,DNele(i)*Nfp,'double');
%            St( offset+1:offset+DhNele(i)*Nfp,3,j)=tmp(1:DhNele(i)*Nfp);
            fclose(fidb);
        end
    end
    offset=offset+DhNele(i)*Nfp;
    offsetV=offsetV+DhNele(i)*pNp;
end
disp('time invariant done')

vtuxml_var(fid,pid);
vtuxml_var(fid,Sni);
vtuxml_var(fid,a  );
for i=1:Nt
    vtuxml_var(fid,Sn(:,1,i)*amp);
end
for i=1:Nt
    vtuxml_var(fid,Vt(:,1,i)*amp);
end
for i=1:Nt
    vtuxml_var(fid,Ut(:,1,i)*amp);
end
for i=1:Nt
    vtuxml_var(fid, f(:,1,i));
end
%for i=1:Nt
%    vtuxml_var(fid, bo(:,1,i));
%end
disp('scalar done')

vtuxml_var(fid,Sti);
for i=1:Nt
    vtuxml_var(fid,Sf(:,:,i));
end
for i=1:Nt
    vtuxml_var(fid,St(:,:,i)*amp);
end
for i=1:Nt
    vtuxml_var(fid,dVt(:,:,i)*amp);
end
disp('vector done')

vtuxml_end(fid);

