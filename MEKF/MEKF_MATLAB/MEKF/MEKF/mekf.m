function [xhist, Phist] = mekf(x0, P0, W, V, rN, whist, yhist, dt)

xhist = zeros(7,size(yhist,2));
xhist(:,1) = x0;

Phist = zeros(6,6,size(yhist,2));
Phist(:,:,1) = P0;

for k = 1:(size(yhist,2)-1)
   
    [xp, A] = prediction(xhist(:,k),whist(:,k), dt);
    
%     %Derivative Check
%     Atest = zeros(6);
%     dx = 1e-6*eye(6);
%     for m = 1:3
%         dxin = [qmult(xhist(1:4,k),[0.5*dx(1:3,m);1]); xhist(5:7,k)];
%         dxp = prediction(dxin,whist(:,k));
%         dqp = qmult(qconj(xp(1:4)),dxp(1:4));
%         Atest(:,m) = [2*dqp(1:3); dxp(5:7)-xp(5:7)]/1e-6;
%     end
%     for m = 4:6
%         dxp = prediction(xhist(:,k)+[0;dx(:,m)],whist(:,k));
%         dqp = qmult(qconj(xp(1:4)),dxp(1:4));
%         Atest(:,m) = [2*dqp(1:3); dxp(5:7)-xp(5:7)]/1e-6;
%     end
    
    Pp = A*Phist(:,:,k)*A' + W;
    
    
    [yp, C] = measurement(xp(1:4),rN);
    
%     %Derivative Check
%     Ctest = zeros(6);
%     dx = 1e-6*eye(6);
%     for m = 1:3
%         dxin = [qmult(xp(1:4),[0.5*dx(1:3,m);1]); xp(5:7)];
%         dyp = measurement(dxin,rN);
%         Ctest(:,m) = (dyp-yp)/1e-6;
%     end
%     for m = 4:6
%         dyp = measurement(xp(:,k)+[0;dx(:,m)],rN);
%         Ctest(:,m) = (dyp-yp)/1e-6;
%     end
    
    %Innovation
%     Qp = qtodcm(xp(1:4));
%     Qm = qtodcm(yhist(1:4,k+1));
%     phiST = vee(logSO3(Qp'*Qm));
    phiST = [2*eye(3), zeros(3,1)]*logq(qmult(qconj(yp(1:4)),yhist(1:4,k+1)));
    z = [phiST; yhist(5:end,k+1) - yp(5:end)];
    S = C*Pp*C' + V;
%     z = [yhist(5:end,k+1) - yp(5:end)];
%     S = C(4:9,:)*Pp*C(4:9,:)' + V(4:9,4:9);

    %Kalman Gain
    L = (Pp*C')/S;
%     L = (Pp*C(4:9,:)')/S;
    
    %Update
    dx = L*z;

    dq = [0.5*dx(1:3); 1-0.125*(dx(1:3)'*dx(1:3))];
    dq = dq/norm(dq);
    
    xhist(1:4,k+1) = qmult(xp(1:4),dq);
    xhist(5:7,k+1) = xp(5:7) + dx(4:6);
    
    Phist(:,:,k+1) = (eye(6)-L*C)*Pp*(eye(6)-L*C)' + L*V*L';
%     Phist(:,:,k+1) = (eye(6)-L*C(4:9,:))*Pp*(eye(6)-L*C(4:9,:))' + L*V(4:9,4:9)*L';

end


end

