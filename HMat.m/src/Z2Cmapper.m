function order = Z2Cmapper(n,minn,tensor)

if max(n) <= minn
    order = reshape(tensor,prod(n),[]);
    return;
end
nl = floor(n/2);
nr = n-nl;
nv = zeros(size(n));
nrange = cell(2,1);
offset = 0;
order = zeros(prod(n),1);
for i = 0:2^2-1
    for d = 0:2-1
        nv(d+1) = nl(d+1)*(1-mod(floor(i/2^d),2))+nr(d+1)*mod(floor(i/2^d),2);
        nrange{d+1} = nl(d+1)*mod(floor(i/2^d),2)+(1:floor(nv(d+1)));
    end
    order(offset+(1:prod(nv))) = Z2Cmapper(nv,minn,tensor(nrange{1},nrange{2}));
    offset = offset + prod(nv);
end

end
