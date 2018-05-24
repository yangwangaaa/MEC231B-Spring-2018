function cost = evaluateWorstCostManyModels(A,E,C,F,Lmat)
m = size(A,3);
cost = inf;
for i = 1:m
    ANew = A(:,:,i);
    ENew = E(:,:,i);
    CNew = C(:,:,i);
    FNew = F(:,:,i);
    cost2 = evaluateCost(ANew,ENew,CNew,FNew,Lmat);
    if cost2 < cost
        cost = cost2;
    else
        % Nothing
    end
end
end