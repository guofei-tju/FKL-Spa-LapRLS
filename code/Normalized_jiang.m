function S=Knormalized_jiang(K)
%kernel normilization
    n=length(K);
    K = K./repmat(sum(K,2),1,n);
    for i=1:n
        K(i,i)=1;
    end
    S=K;
end