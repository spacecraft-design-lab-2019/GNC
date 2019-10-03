## note: the quaternion in all of these functions goes:
## [scalar; vector]

function quaternion_euler(q)
    #this function converts a quaternion into its euler angles

    euler = [atan(2*(q[1]*q[2] + q[3]*q[4]),1-2*(q[2]^2+q[3]^2));
            asin(2*(q[1]*q[3] - q[4]*q[2]));
            atan(2*(q[1]*q[4] + q[2]*q[3]),1-2*(q[3]^2+q[4]^2))]

    return euler

end

function euler_quaternion(euler)
    #this function (poorly) converts a set of euler angles into a quaternion

    ϕ = euler[1]
    θ = euler[2]
    ψ = euler[3]

    quat = [cos(ϕ/2)*cos(θ/2)*cos(ψ/2) +
    sin(ϕ/2)*sin(θ/2)*sin(ψ/2);
    sin(ϕ/2)*cos(θ/2)*cos(ψ/2) -
    cos(ϕ/2)*sin(θ/2)*sin(ψ/2);
    cos(ϕ/2)*sin(θ/2)*cos(ψ/2) +
    sin(ϕ/2)*cos(θ/2)*sin(ψ/2);
    cos(ϕ/2)*cos(θ/2)*sin(ψ/2) -
    sin(ϕ/2)*sin(θ/2)*cos(ψ/2)]

    return quat

end

function quaternion_rot(q)
    # this function takes in a quaternion and spits out a rotation matrix
    w = q[1]
    x = q[2]
    y = q[3]
    z = q[4]
    ROT = [1-2y^2-2z^2 2x*y-2z*w 2x*z+2y*w;
       2x*y+2z*w 1-2x^2-2z^2 2y*z-2x*w;
       2x*z-2y*w 2y*z+2x*w 1-2x^2-2y^2]

end

function rot_quaternion(ROT)
    # this function takes in a rotation matrix and spits out a quaternion

    # w= sqrt(1 + ROT[1,1] + ROT[2,2] + ROT[3,3])/2
    # x = (ROT[3,2] - ROT[2,3])/(4*w)
    # y = (ROT[1,3] - ROT[3,1])/(4*w)
    # z = (ROT[2,1] - ROT[1,2])/(4*w)
    #
    # q = [w;x;y;z]

    r11 = ROT[1,1]
    r12 = ROT[1,2]
    r13 = ROT[1,3]
    r21 = ROT[2,1]
    r22 = ROT[2,2]
    r23 = ROT[2,3]
    r31 = ROT[3,1]
    r32 = ROT[3,2]
    r33 = ROT[3,3]


    # shepperd's method

    i = sortperm([(r11+r22+r33);r11;r22;r33])[end]

    if i == 1
        r1 = sqrt(1+r11+r22+r33)
        q = .5*[r1;(r32-r23)/(r1);(r13-r31)/r1;(r21-r12)/r1]
    elseif i == 2
        r2 = sqrt(1+r11-r22-r33)
        q = .5*[(r32-r23)/(r2);r2;(r12+r21)/r2;(r31+r13)/r2]
    elseif i == 3
        r3 = sqrt(1-r11+r22-r33)
        q = .5*[(r13-r31)/r3;(r12+r21)/r3;r3;(r23+r32)/r3]
    elseif i == 4
        r4 = sqrt(1-r11-r22+r33)
        q = .5*[(r21-r12)/r4;(r31+r13)/r4;(r32+r23)/r4;r4]
    else
        println("error in Shepperd's method")
    end

    return q

end
