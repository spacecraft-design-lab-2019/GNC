clear all

load mekf_inputs

rN = [rN1, rN2];
%yhist = [qSThist; rB1hist; rB2hist];
yhist = [rB1hist; rB2hist];

q0 = triad(rN1,rN2,rB1hist(:,1),rB2hist(:,2));
%q0 = yhist(1:4,1);

x0 = [q0; 0; 0; 0];
P0 = blkdiag((.1*pi/180)^2*eye(3), (1*pi/180)^2*eye(3)); %10 deg. and 10 deg/sec 1-sigma uncertainty

[xhist,Phist] = mekf(x0,P0,W,V,rN,whist,yhist,dt);

load mekf_truth

e = zeros(size(qtrue));
for k = 1:size(qtrue,2)
    e(:,k) = qmult(qconj(qtrue(:,k)), xhist(1:4,k));
end

close all

figure(1);
subplot(4,1,1);
plot(qtrue(1,:));
hold on;
plot(xhist(1,:));
title('Attitude');
legend('True', 'Estimated');
subplot(4,1,2);
plot(qtrue(2,:));
hold on;
plot(xhist(2,:));
subplot(4,1,3);
plot(qtrue(3,:));
hold on;
plot(xhist(3,:));
subplot(4,1,4);
plot(qtrue(4,:));
hold on;
plot(xhist(4,:));

figure(2);
subplot(3,1,1);
plot((360/pi)*e(1,:));
hold on
plot((360/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
plot(-(360/pi)*sqrt(squeeze(Phist(1,1,:))),'r');
title('Attitude Error');
subplot(3,1,2);
plot((360/pi)*e(2,:));
hold on
plot((360/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
plot(-(360/pi)*sqrt(squeeze(Phist(2,2,:))),'r');
ylabel('degrees');
subplot(3,1,3);
plot((360/pi)*e(3,:));
hold on
plot((360/pi)*sqrt(squeeze(Phist(3,3,:))),'r');
plot(-(360/pi)*sqrt(squeeze(Phist(3,3,:))),'r');

figure(3);
subplot(3,1,1);
plot(xhist(5,:)-btrue(1,:));
hold on
plot(2*sqrt(squeeze(Phist(4,4,:))),'r');
plot(-2*sqrt(squeeze(Phist(4,4,:))),'r');
title('Bias Error');
subplot(3,1,2);
plot(xhist(6,:)-btrue(2,:));
hold on
plot(2*sqrt(squeeze(Phist(5,5,:))),'r');
plot(-2*sqrt(squeeze(Phist(5,5,:))),'r');
subplot(3,1,3);
plot(xhist(7,:)-btrue(3,:));
hold on
plot(2*sqrt(squeeze(Phist(6,6,:))),'r');
plot(-2*sqrt(squeeze(Phist(6,6,:))),'r');

