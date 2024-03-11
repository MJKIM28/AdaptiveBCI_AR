function D2 = distfun(ZI,ZJ)
for j = 1:size(ZJ,1)
y = ZJ(j,:);
% D = mysigmoid(ZI,ZI) + mysigmoid(y,y) - 2*mysigmoid(ZI,y);

D = 1 - mysigmoid(ZI,y)/(sqrt(mysigmoid(ZI,ZI)*mysigmoid(y,y)));

D2(j,1) = D;
end

% ZInew = repmat(ZI,size(ZJ,1),1);
% Dnew = mysigmoid(ZInew,ZInew) + mysigmoid(ZJ,ZJ) - 2*mysigmoid(ZInew,ZJ);
% D2 = diag(Dnew);

end