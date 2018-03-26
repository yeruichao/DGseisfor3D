function snap2vtu3d_new(configfile,v0,basename)

[cond,scalvar,vectvar,snapname,outfile]=read_visual_para(configfile);

fname=[snapname,'_Para.txt'];
Para=load(fname);
Nbele =Para(1); % glob_Nele
Np    =Para(2); % pNp
Nsele =Para(4); % Nsele
Nproc =Para(5); % Nprocs
Nnod  = Nbele*Np;
Nele  = Nbele*Nsele;
DNele = zeros(Nproc,1);
DNNele = zeros(Nproc,1);

Nscalvar=size(scalvar,1);Nscal=0;
Nvectvar=size(vectvar,1);Nvect=0;

for i=1:Nscalvar
    Nscal=Nscal+scalvar{i}.final-scalvar{i}.start+1;
end
for i=1:Nvectvar
    Nvect=Nvect+vectvar{i}.final-vectvar{i}.start+1;
end
Scalname=cell(Nscal,1);Vectname=cell(Nvect,1);

iscal=0;
for i=1:Nscalvar
    if(scalvar{i}.start < 0)
        iscal=iscal+1;
        Scalname{iscal}=scalvar{i}.varname;
    else
        for j=scalvar{i}.start:scalvar{i}.final
            iscal=iscal+1;
            Scalname{iscal}=[scalvar{i}.varname,'_t',num2str(j,'%d')];
        end
    end
end

ivect=0;
for i=1:Nvectvar
    if(vectvar{i}.start < 0)
        ivect=ivect+1;
        Vectname{ivect}=vectvar{i}.varname{1}(1:end-1);
    else
        for j=vectvar{i}.start:vectvar{i}.final
            ivect=ivect+1;
            Vectname{ivect}=[vectvar{i}.varname{1}(1:end-1),...
                '_t',num2str(j,'%d')];
        end
    end
end

if(nargin>=2)
    globID=zeros(Nbele,1);
end

offset=0;
for i=1:Nproc
    fname=[snapname,'_Subdomain_',num2str(i-1),'.dat'];
    fidb=fopen(fname,'r');
    DNNele(i)=fread(fidb,1,'int');            % local_Nele
    DNele(i)=fread(fidb,1,'int');   % local_Nhele
    if(nargin>=2)
        globID(offset+1:offset+DNele(i))=fread(fidb,DNele(i),'int'); % globID
        offset=offset+DNele(i);
    end
    fclose(fidb);
end
if(sum(DNele)~=Nbele)
    disp('Something wrong with the number of elements.');
end

p=zeros(Nnod,3);offset=0;
for i=1:Nproc
    fname=[snapname,'_X_',num2str(i-1),'.dat'];
    fidb=fopen(fname,'r');
    p(offset+1:offset+DNele(i)*Np,1)=fread(fidb,DNele(i)*Np,'double');
    fclose(fidb);
    fname=[snapname,'_Y_',num2str(i-1),'.dat'];
    fidb=fopen(fname,'r');
    p(offset+1:offset+DNele(i)*Np,2)=fread(fidb,DNele(i)*Np,'double');
    fclose(fidb);
    fname=[snapname,'_Z_',num2str(i-1),'.dat'];
    fidb=fopen(fname,'r');
    p(offset+1:offset+DNele(i)*Np,3)=fread(fidb,DNele(i)*Np,'double');
    fclose(fidb);
    offset=offset+DNele(i)*Np;
end

disp('Coordinate done')

t=zeros(4*Nele,1);offset=0;counter=0;
for i=1:Nproc
    fname=[snapname,'_Conn_',num2str(i-1),'.dat'];
    fidb=fopen(fname,'r');
    tmp=fread(fidb,DNNele(i)*Nsele*4,'int')+counter;
    tmp=reshape(tmp,DNNele(i)*Nsele,4);tmp=tmp(1:DNele(i)*Nsele,:)';
    t(offset+1:offset+DNele(i)*Nsele*4)=tmp(:);
    offset=offset+DNele(i)*Nsele*4;counter=counter+DNele(i)*Np;
    fclose(fidb);
end
t=reshape(t,4,Nele);

disp('Connection done')

if(nargin>=2)
    globID=ones(Nsele,1)*globID';
    globID=globID(:);
    if(not(isempty(v0)))
        v0=v0(globID);
    else
        v0=[];
    end
end

%[t,nid]=pickeles(p,t,cond,v0,Nsele);
%[t,nid]=pickeles(p,t,cond,v0);
[t,nid]=pickeles_new(p,t,cond,v0,Nsele,basename,Np,globID);

t=t'-1;
p=p(nid,:);
nnode=size(p,1);
nnele=size(t,2);

fid=vtuxml_head(outfile,nnode,nnele,Scalname,Vectname);
vtuxml_xyz(fid,p);
vtuxml_conn(fid,t);

clear t p;

disp(['Original number of nodes: ',num2str(Nnod,'%d')]);
disp(['Sampled number of nodes: ',num2str(nnode,'%d'), ...
    ', ',num2str(nnode/Nnod*100,'%f'),' percent of total nodes']);
disp(['Original number of elements: ',num2str(Nele,'%d')]);
disp(['Sampled number of elements: ',num2str(nnele,'%d'), ...
    ', ',num2str(nnele/Nele*100,'%f'),' percent of total elements']);

for i=1:Nscalvar
    if(scalvar{i}.start < 0)
        cname=[snapname,'_',scalvar{i}.varname];
        p=zeros(Nnod,1);offset=0;
        for j=1:Nproc
            fname=[cname,'_',num2str(j-1,'%d'),'.dat'];
            fidb=fopen(fname,'r');
            p(offset+1:offset+DNele(j)*Np)= ...
                fread(fidb,DNele(j)*Np,'double');
            offset=offset+DNele(j)*Np;
            fclose(fidb);
        end
        disp(['Value of ',scalvar{i}.varname,': ',num2str(min(p)), ...
            ' ~ ',num2str(max(p))]);
        p=p(nid);
        vtuxml_var(fid,p);
        clear p;
        disp(['Writing ',scalvar{i}.varname,' done.'])
    else
        for isnap=scalvar{i}.start:scalvar{i}.final
            cname=[snapname,'_',scalvar{i}.varname,'_',num2str(isnap,'%d')];
            p=zeros(Nnod,1);offset=0;
            for j=1:Nproc
                fname=[cname,'_',num2str(j-1,'%d'),'.dat'];
                fidb=fopen(fname,'r');
                p(offset+1:offset+DNele(j)*Np)= ...
                    fread(fidb,DNele(j)*Np,'double');
                offset=offset+DNele(j)*Np;
                fclose(fidb);
            end
            disp(['Value of ',scalvar{i}.varname,' snapshot No. ',...
                num2str(isnap,'%d'),': ',num2str(min(p)), ...
                ' ~ ',num2str(max(p))]);
            p=p(nid);
            vtuxml_var(fid,p);
            clear p;
            disp(['Writing ',scalvar{i}.varname,' snapshot No. ', ... 
                num2str(isnap,'%d'), ' done.']);
        end
    end
end

for i=1:Nvectvar
    if(vectvar{i}.start < 0)
        v=zeros(Nnod,3);
        for k=1:3
            offset=0;
            cname=[snapname,'_',vectvar{i}.varname{k}];
            for j=1:Nproc
                fname=[cname,'_',num2str(j-1),'.dat'];
                fidb=fopen(fname,'r');
                v(offset+1:offset+DNele(j)*Np,k)= ...
                    fread(fidb,DNele(j)*Np,'double');
                fclose(fidb);
                offset=offset+DNele(j)*Np;
            end
        end
        disp(['Value of ',vectvar{i}.varname{i}(1:end-1),': ', ...
            num2str(min(v)),' ~ ',num2str(max(v))]);
        v=v(nid,:);
        vtuxml_var(fid,v);
        clear v;
        disp(['Writing ',vectvar{i}.varname{i}(1:end-1),' done.'])
    else
        for isnap=vectvar{i}.start:vectvar{i}.final
            v=zeros(Nnod,3);
            for k=1:3
                offset=0;
                cname=[snapname,'_',vectvar{i}.varname{k},'_', ...
                    num2str(isnap,'%d')];
                for j=1:Nproc
                    fname=[cname,'_',num2str(j-1),'.dat'];
                    fidb=fopen(fname,'r');
                    v(offset+1:offset+DNele(j)*Np,k)= ...
                        fread(fidb,DNele(j)*Np,'double');
                    fclose(fidb);
                    offset=offset+DNele(j)*Np;
                end
            end
            disp(['Value of ',vectvar{i}.varname{i}(1:end-1),...
                ' snapshot No. ',num2str(isnap,'%d'),' : ', ...
                num2str(min(v)),' ~ ',num2str(max(v))]);
            v=v(nid,:);
            vtuxml_var(fid,v);
            clear v;
            disp(['Writing ',vectvar{i}.varname{i}(1:end-1), ...
                ' snapshot No. ',num2str(isnap,'%d'),' done.'])
        end
    end
end

vtuxml_end(fid);

end
