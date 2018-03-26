function tt=trisubele(N)

    if N==1
	tt=1:3;
	return
    end
    Nsele = N*(N+1)/2 + (N-1)*N/2;
    tt=zeros(Nsele,3);
    M=zeros(N+1,N+1);

    sk=0;
    for k=1:N+1
        for j=1:N+2-k
            sk=sk+1;
            M(j,k)=sk;
        end
    end
    sk=0;
    for k=1:N
        for j=1:N+1-k
            sk=sk+1;
            tt(sk,1)=M(j,k)  ;
            tt(sk,2)=M(j+1,k);
            tt(sk,3)=M(j,k+1);
        end
    end
    for k=1:N-1
        for j=1:N-k
            sk=sk+1;
            tt(sk,1)=M(j+1,k)  ;
            tt(sk,2)=M(j+1,k+1);
            tt(sk,3)=M(j,k+1)  ;
        end
    end

end

