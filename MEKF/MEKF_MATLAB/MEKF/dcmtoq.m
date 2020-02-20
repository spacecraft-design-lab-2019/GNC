function q = dcmtoq(R)

    if ndims(R) == 2
        q = Rtoq(R);
    else %ndims(R) == 3
        q = zeros(4,size(R,3));
        q(:,1) = Rtoq(R(:,:,1));
        for k = 2:size(R, 3)
            q(:,k) = Rtoq(R(:,:,k));
            if norm(q(:,k) - q(:,k-1)) > norm(q(:,k) + q(:,k-1))
                q(:,k) = -q(:,k);
            end
        end
    end
    
    function q = Rtoq(R)
        q = [0 0 0 0]';
        q(4) = .5*sqrt(1+trace(R));
        q(1) = sign(R(3,2)-R(2,3))*abs(.5*sqrt(1+R(1,1)-R(2,2)-R(3,3)));
        q(2) = sign(R(1,3)-R(3,1))*abs(.5*sqrt(1-R(1,1)+R(2,2)-R(3,3)));
        q(3) = sign(R(2,1)-R(1,2))*abs(.5*sqrt(1-R(1,1)-R(2,2)+R(3,3)));
        
        %Ensure normalization
        q = q/sqrt(q'*q);
    end

end