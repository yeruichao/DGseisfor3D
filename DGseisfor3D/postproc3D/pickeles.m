function IDm=pickeles(p,cond)

pNp=size(p,1);
Nele=size(p,2);

IDm=false(Nele,1);
Ncond=length(cond);
for i=1:Ncond
    ID=true(Nele,1);
    Nscond=length(cond{i});
    for j=1:Nscond
        if(abs(cond{i}{j}.type) == 1)
            n=cond{i}{j}.value;
            n=n(:)/norm(n);
            p0n=reshape(reshape(p,pNp*Nele,3)*n,pNp,Nele);
            if(cond{i}{j}.type < 0)
                ID=ID & any(p0n <= cond{i}{j}.offset,1)';
            else
                ID=ID & any(p0n >= cond{i}{j}.offset,1)';
            end
        elseif(cond{i}{j}.type == 0)
	    a=(6378.136/6371)^-2;
	    b=(6356.752/6371)^-2;
            r=sqrt(a*(p(:,:,1)-cond{i}{j}.value(1)).^2 ...
                  +a*(p(:,:,2)-cond{i}{j}.value(2)).^2 ...
                  +b*(p(:,:,3)-cond{i}{j}.value(3)).^2);
            if(cond{i}{j}.offset < 0)
                ID=ID & any(r <= -cond{i}{j}.offset,1)';
            else
                ID=ID & any(r >= cond{i}{j}.offset,1)';
            end
        end
    end
    IDm=IDm | ID;
end

end
