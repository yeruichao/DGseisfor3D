function P=showsnap2D(basename,varname)

filename=strcat('../snapshot2/',basename,'_Para.txt');
disp(filename)
Para=load(filename);
n1=Para(1);
n2=Para(2);
filename=strcat('../snapshot2/',basename,'_',varname);
fid=fopen(filename,'r');
P=fread(fid,n1*n2,'float');
P=reshape(P,n1,n2);
imagesc(P);
fclose(fid);
end