function snap2vtu3d(configfile,att)

[cond,regions,vars,snapname,outfile,tmpdir]=read_visual_para(configfile);

fname=[snapname,'_Para.txt'];
if(not(exist(fname,'file')))
    disp(fname)
end
Para=load(fname);
glob_Nele = Para(1); % glob_Nele
pNp   = Para(2); 
Nfp   = Para(3); 
Nproc = Para(4); 
Np    = pNp*3/Nfp-3

disp(['glob_Nele=',num2str(glob_Nele,'%d')])

Nsnap = length(vars.snap);

Nscal=0;
for i=1:size(vars.scalvar,1)
    if(vars.scalvar{i}.flag)
        Nscal=Nscal+1;
    end
end
for i=1:size(vars.scalvarT,1)
    if(vars.scalvarT{i}.flag)
        Nscal=Nscal+Nsnap;
    end
end

Nvect=0;
for i=1:size(vars.vectvar,1)
    if(vars.vectvar{i}.flag)
        Nvect=Nvect+1;
    end
end
for i=1:size(vars.vectvarT,1)
    if(vars.vectvarT{i}.flag)
        Nvect=Nvect+Nsnap;
    end
end

scalname=cell(Nscal,1);vectname=cell(Nvect,1);

iscal=0;
for i=1:size(vars.scalvar,1)
    if(vars.scalvar{i}.flag)
        iscal=iscal+1;
        scalname{iscal}=vars.scalvar{i}.varname;
    end
end
for i=1:size(vars.scalvarT,1)
    if(vars.scalvarT{i}.flag)
        for j=1:Nsnap
            iscal=iscal+1;
            scalname{iscal}=[vars.scalvarT{i}.varname,...
                '_t',num2str(vars.snap(j),'%04d')];
        end
    end
end

ivect=0;
for i=1:size(vars.vectvarT,1)
    if(vars.vectvarT{i}.flag)
        for j=1:Nsnap
            ivect=ivect+1;
            vectname{ivect}=[vars.vectvarT{i}.varname,...
                '_t',num2str(vars.snap(j),'%04d')];
        end
    end
end
for i=1:size(vars.vectvar,1)
    if(vars.vectvar{i}.flag)
        ivect=ivect+1;
        vectname{ivect}=vars.vectvar{i}.varname;
    end
end

DNele  = zeros(Nproc,1);
Dsele  = zeros(Nproc,1);
Dmask  =  cell(Nproc,1);
Dcoord =  cell(Nproc,1);
Dpml   = false(Nproc,1);
k_reg = not(isempty(regions)) && Nargin>=2 && not(isempty(att));

for i=1:Nproc
    fname=[snapname,'_Subdomain_',num2str(i-1),'.dat'];
    if(not(exist(fname,'file')))
        disp([fname,' does not exist']);return
    end
    fidb=fopen(fname,'r');
    Nele   = fread(fidb,1,'int');  % local_Nele
    Nhele  = fread(fidb,1,'int');  % local_Nhele
    Ninele = fread(fidb,1,'int');  % local_Ninele
    kpml   = fread(fidb,1,'int');  % pml flag
    xyz = fread(fidb,pNp*Nele*3,'double');  % nodal coordinates
    xyz = reshape(xyz,pNp,Nele,3);
    if(k_reg)
        mask = false(Nele,1);
        globID=fread(fidb,Nele,'int'); % globID
        loc_att=att(globID);
        for j=1,length(regions)
            mask = mask | (loc_att==regions(j));
        end
    else
        mask = true(Nele,1);
    end
    if(Nele>Nhele)
        mask(Nhele+1:end)=false;
    end
    fclose(fidb);
    mask = mask & pickeles(xyz,cond);
    Dmask{i}=find(mask);
    DNele(i)=Nele;
    Dpml(i)=kpml~=0;
    Dsele(i)=length(Dmask{i});
    if(Dsele(i)>0)
        Dcoord{i}=xyz(:,mask,:);
    end
end

tt=tetsubele(Np);
Nsele=size(tt,1);
tt=tt';

Nele=sum(Dsele);
vis_Nnod=Nele*pNp;
vis_Nele=Nele*Nsele;

fid=vtuxml_head(outfile,vis_Nnod,vis_Nele,scalname,vectname);
fclose(fid);

p=zeros(pNp,Nele,3);
offset=0;
for i=1:Nproc
    if(Dsele(i)<=0)
        continue;
    end
    p(:,offset+1:offset+Dsele(i),:)=Dcoord{i};
    offset=offset+Dsele(i);
end
p=reshape(p,vis_Nnod,3);

fid=fopen([outfile,'.vtu'],'a');
vtuxml_xyz(fid,p);
disp('coordinate done')
clear p;

t=pNp*ones(4*Nsele,1)*(0:Nele-1);
t=t+tt(:)*ones(1,Nele);
t=reshape(t,4,vis_Nele)-1;

vtuxml_conn(fid,t);
disp('connection done')
clear t;

disp(['Original number of elements: ',num2str(glob_Nele*Nsele,'%d')]);
disp(['Sampled number of elements: ',num2str(vis_Nele,'%d'), ...
    ', (',num2str(vis_Nele/glob_Nele/Nsele*100,'%f'),' percent).']);

% Cp
if(vars.scalvar{1}.flag)
    p=read_sample_file(snapname,'Cp',DNele,Dsele,Dmask,...
        pNp,Nproc,1,0,[],true(Nproc,1));
    p=reshape(p,vis_Nnod,1);
    vtuxml_var(fid,p);
end

% Cs
if(vars.scalvar{2}.flag)
    p=read_sample_file(snapname,'Cs',DNele,Dsele,Dmask,...
        pNp,Nproc,1,0,[],true(Nproc,1));
    p=reshape(p,vis_Nnod,1);
    vtuxml_var(fid,p);
end

% rho
if(vars.scalvar{3}.flag)
    p=read_sample_file(snapname,'Material',DNele,Dsele,Dmask,...
        pNp,Nproc,1,1,'int',true(Nproc,1));
    p=reshape(p,vis_Nnod,1);
    vtuxml_var(fid,p);
end

% wavefield: V1,V2,V3,E11,E22,E33,E23,E13,E12
for i=1:Nsnap
    isnap=vars.snap(i);
    varname=['Wave_',num2str(isnap,'%d')];
    p=read_sample_file(snapname,varname,DNele,Dsele,Dmask,...
        pNp,Nproc,9,0,[],true(Nproc,1));
    tmpfid=fopen([tmpdir,varname,'.V.tmpfile'],'w');
    fwrite(tmpfid,p(:,:,1:3),'double');
    fclose(tmpfid);
    for j=1:6
        tmpfid=fopen([tmpdir,varname,'.E',num2str(j,'%d'),'.tmpfile'],'w');
        fwrite(tmpfid,p(:,:,j+3),'double');
        fclose(tmpfid);
    end
end

% E
for i=1:6
    if(vars.scalvarT{i}.flag)
        for j=1:Nsnap
            isnap=vars.snap(j);
            varname=['Wave_',num2str(isnap,'%d')];
            tmpfid=fopen([tmpdir,varname,'.E',num2str(i,'%d'),'.tmpfile'],'r');
            p=fread(tmpfid,vis_Nnod,'double');
            vtuxml_var(fid,p);
        end
    end
end

% V
if(vars.vectvarT{1}.flag)
    for j=1:Nsnap
        isnap=vars.snap(j);
        varname=['Wave_',num2str(isnap,'%d')];
        tmpfid=fopen([tmpdir,varname,'.V.tmpfile'],'r');
        p=fread(tmpfid,vis_Nnod*3,'double');
        p=reshape(p,vis_Nnod,3);
        vtuxml_var(fid,p);
    end
end

unix(['rm   ',tmpdir,'*.tmpfile']);

% PML
if(vars.vectvar{1}.flag)
    p=read_sample_file(snapname,'PMLdamp',DNele,Dsele,Dmask,...
        pNp,Nproc,3,0,[],Dpml);
    p=reshape(p,vis_Nnod,3);
    vtuxml_var(fid,p);
end

vtuxml_end(fid);

end


function p=read_sample_file(snapname,varname,DNele,Dsele,Dmask,...
    pNp,Nproc,Ncomp,offset,offset_type,iDomain)

    Nele=sum(Dsele);
    p=zeros(pNp,Nele,Ncomp);
    boffset=0;
    for i=1:Nproc
        if(Dsele(i)<=0)
            continue;
        end
        fname=[snapname,'_',varname,'_',num2str(i-1),'.dat'];
        if(not(iDomain(i)) || not(exist(fname,'file')))
            disp([fname,' does not exist']);
            tmp=zeros(pNp,DNele(i),Ncomp);
            p(:,boffset+1:boffset+Dsele(i),:)=tmp(:,Dmask{i},:);
        else
            fidb=fopen(fname,'r');
            if(offset>0)
                fread(fidb,offset,offset_type);
            end
            tmp=fread(fidb,pNp*DNele(i)*Ncomp,'double');
            tmp=reshape(tmp,pNp,DNele(i),Ncomp);
            p(:,boffset+1:boffset+Dsele(i),:)=tmp(:,Dmask{i},:);
        end
        boffset=boffset+Dsele(i);
    end
end
