function [p,t,pt,at,neigh]=read_mesh3d(fname)

    fname1=fname;
    fid = fopen([fname1 '.ele'],'r'); % read in .ele file
    A=fscanf(fid,'%d %d %d',3);
    Nele = A(1); % Number of tetrahedron
    Nattr = A(3); % Number of attributes
    A = fscanf(fid,'%d',Nele*(5+Nattr));
    A = reshape(A,5+Nattr,Nele)';
    t=A(:,[2 3 4 5]); % tri ID, node 1~3
    fclose(fid);
    if(Nattr>0)
        at=A(:,6);
    else
        at=[];
    end
    fid = fopen([fname1 '.node'],'r'); %read in .node file
    A=fscanf(fid,'%d',4); 
    Nnodes = A(1); % Number of nodes
    Nattr = A(3); % Number of attributes
    Nbdry = A(4); % Number of boundaries
    A=fscanf(fid,'%g',Nnodes*(4+Nattr+Nbdry));
    A = reshape(A,4+Nattr+Nbdry,Nnodes)';
    p = A(:,[2 3 4]);
    pt=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:)+p(t(:,4),:))/4;
    
    fname1=[fname,'.neigh'];
    if(exist(fname1,'file'))
        fid=fopen(fname1,'r'); %read in .neigh file
        fscanf(fid,'%d',2);
        A = fscanf(fid,'%d',Nele*5);
        A = reshape(A,5,Nele)';
        neigh=A(:,2:5);
    else
        neigh=[];
    end
end

