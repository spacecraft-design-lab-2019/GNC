function [rN1, rN2, rB1hist, rB2hist, qSThist, rB1true, rB2true] = sensors(xhist,V1,V2,V3)

sqrtV1 = chol(V1);
sqrtV2 = chol(V2);
sqrtV3 = chol(V3);

qhist = xhist(1:4,:);

rN1 = randn(3,1);
rN1 = rN1/norm(rN1);

rN2 = randn(3,1);
rN2 = rN2/norm(rN2);

rB1true = zeros(3,size(qhist,2));
rB2true = zeros(3,size(qhist,2));

rB1hist = zeros(3,size(qhist,2));
rB2hist = zeros(3,size(qhist,2));

qSThist = zeros(size(qhist));

for k = 1:size(qhist,2)
    QBN = qtodcm(qhist(:,k))';
    rB1true(:,k) = QBN*rN1;
    rB2true(:,k) = QBN*rN2;
    
    %Add noise and re-normalize
    rB1hist(:,k) = rB1true(:,k) + sqrtV1*randn(3,1);
    rB1hist(:,k) = rB1hist(:,k)/norm(rB1hist(:,k));
    rB2hist(:,k) = rB2true(:,k) + sqrtV2*randn(3,1);
    rB2hist(:,k) = rB2hist(:,k)/norm(rB2hist(:,k));
    
    phi = sqrtV3*randn(3,1);
    dq = [phi/2; 1-.125*(phi'*phi)];
    dq = dq/norm(dq);
    qSThist(:,k) = qmult(qhist(:,k),dq);
end


end

