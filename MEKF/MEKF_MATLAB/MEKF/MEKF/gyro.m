function [whist,bhist] = gyro(xhist,b0,Vw,Vb)

sqrtVw = chol(Vw);
sqrtVb = chol(Vb);

wtrue = xhist(5:7,:);

bhist = zeros(3,size(wtrue,2));
bhist(:,1) = b0;

whist = zeros(3,size(wtrue,2));
whist(:,1) = wtrue(:,1) + bhist(:,1) + sqrtVw*randn(3,1);

for k = 2:size(xhist,2)
    bhist(:,k) = bhist(:,k-1) + sqrtVb*randn(3,1);
    whist(:,k) = wtrue(:,k) + bhist(:,k) + sqrtVw*randn(3,1);
end

end

