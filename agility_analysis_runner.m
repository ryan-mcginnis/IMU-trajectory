%% Agility Run Analysis Runner
% This script provides implementation of code for estimating velocty and 
% displacement from IMU data when trajectory has known waypoints.
% Written by Ryan S. McGinnis - ryan.mcginis14@gmail.com - July 9, 2016

% Copyright (C) 2016  Ryan S. McGinnis
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% Method used in:
% McGinnis, Ryan S., et al. "Inertial sensor and cluster analysis for 
% discriminating agility run technique." IFAC-PapersOnLine 48.20 (2015): 
% 423-428.


%% Load data file and import orientation functions
addpath('C:\Users\ryanmcg\Documents\Repos\IMU-orientation'); % add location of IMU orientaiton repo to path
load('example_data.mat');
t = data.t;
a = data.a; %m/s^2
w = data.w; %rad/s
spine_dir = data.spine_dir;
ind = data.ind_still; % indices for still section at beginning of trial

% Manage memory
clear data


%% Extract raw kinematics
%Define orientation 
q = get_orientation_compfilter_quaternion(t, a, w * 180/pi, ind);

%Define variables for calculation
t_j = t;
a_j = quaternRot(q, a) - (t_j.^0)*[0,0,norm(mean(a(ind,:)))];
w_j = quaternRot(q, w);
vert = quaternRot(q, spine_dir); %orientation of spine in lab frame
ant = quaternRot(q,[0,0,-1]); %orientation of anterior direction in lab frame - assumes known posterior device placement 

%Resolve acceleration into a-p and m-l components
%anterior-posterior accel
ant_h = [ant(:,1)./sqrt(sum(ant(:,1:2).^2,2)), ant(:,2)./sqrt(sum(ant(:,1:2).^2,2))];
a_apT = dot(a_j(:,1:2),ant_h,2);

%medio-lateral accel
ml_h = cross(ant,(ant(:,1).^0)*[0,0,1]);
a_mlT = dot(a_j(:,1:2),ml_h(:,1:2),2);

%Define tilt, rom, and azimuthal angle
%tilt
tilt = real(acosd(dot(vert,(vert(:,1).^0)*[0,0,1],2)));

%rom
rom = max(tilt)-min(tilt);

%azimuthal angle
aang = ang_scan(real(atan2d(ant_h(:,1),ant_h(:,2)))); 
aang = aang - aang(1);

%Resolve tilt into AP and ML components
th = vert - dot(vert,(vert(:,1).^0)*[0,0,1],2)*[0,0,1];
tap = real(atand(dot(th(:,1:2),ant_h,2)./dot(vert,(vert(:,1).^0)*[0,0,1],2)));
tml = real(atand(dot(th,ml_h,2)./dot(vert,(vert(:,1).^0)*[0,0,1],2)));

%Filter out noise in data due to foot falls
fs = 1/mean(diff(t_j)); % sampling frequency
filter_order = 2;   % 2nd order filter
cut_off = 1.5*(5/(t_j(end)-t_j(1))); % cutoff frequency (Hz) - based on number of turns in course and total time to complete
[num,den] = butter(filter_order,cut_off/(fs/2),'low');

a_mlF = filtfilt(num,den,a_mlT);
t_mlF = filtfilt(num,den,tml);
a_apF = filtfilt(num,den,a_apT);
t_apF = filtfilt(num,den,tap);
w_zF = filtfilt(num,den,w_j(:,3));


%% Find turns based on filtered ml acceleration
%Center ml acceleration
c = [t_j.^0, t_j] \ a_mlF;
fit_c = [t_j.^0, t_j]*c;
a_mlFc = a_mlF-fit_c;

%Find points spanning peaks
p = ones(5,2);
ms = ones(5,1);
thresh = 0.9*rms(a_mlFc);
p(1,1) = find(abs(a_mlFc)>thresh,1,'first');
p(1,2) = find(abs(a_mlFc)<thresh&t_j>t_j(p(1,1)),1,'first');
p(2,1) = find(abs(a_mlFc)>thresh&t_j>t_j(p(1,2)),1,'first');
p(2,2) = find(abs(a_mlFc)<thresh&t_j>t_j(p(2,1)),1,'first');
p(3,1) = find(abs(a_mlFc)>thresh&t_j>t_j(p(2,2)),1,'first');
p(3,2) = find(abs(a_mlFc)<thresh&t_j>t_j(p(3,1)),1,'first');
p(4,1) = find(abs(a_mlFc)>thresh&t_j>t_j(p(3,2)),1,'first');
p(4,2) = find(abs(a_mlFc)<thresh&t_j>t_j(p(4,1)),1,'first');
p(5,1) = find(abs(a_mlFc)>thresh&t_j>t_j(p(4,2)),1,'first');
try
    p(5,2) = find(abs(a_mlFc)<thresh&t_j>t_j(p(5,1)),1,'first');
catch
    p(5,2) = length(t_j);
end
for k=1:5
    [~,ms(k)]=max(abs(a_mlFc(p(k,1):p(k,2))));
    ms(k) = ms(k)+p(k,1)-1;
end
ind_turn = ms;


%% Plot kinematics with identified turns
figure;
set(gcf,'name','ML acceleration');
hold on;
plot(t_j,[t_j.^0, t_j]*c,'k');
plot(t_j,a_mlT,'color',[0.6,0.6,0.6]);
plot(t_j(ind_turn),a_mlF(ind_turn),'ok','MarkerFaceColor','k')
plot(t_j,a_mlF);
plot(t_j,[t_j.^0, t_j]*c + thresh,'--k');
plot(t_j,[t_j.^0, t_j]*c - thresh,'--k');
xlabel('Time (s)','fontsize',16);
ylabel('M-L Accel (m/s^2)','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','ML tilt');
hold on;
plot(t_j,tml,'color',[0.6,0.6,0.6]);
plot(t_j(ind_turn),t_mlF(ind_turn),'.k')
plot(t_j,t_mlF);

figure;
set(gcf,'name','Vertical Angular Velocity');
hold on;
plot(t_j,w_j(:,3),'color',[0.6,0.6,0.6]);
plot(t_j(ind_turn),w_zF(ind_turn),'.k')
plot(t_j,w_zF);


%% Calculate uncorrected velocity and displacement
%Filter out noise in data due to foot falls
fs = 1/mean(diff(t_j)); % sampling frequency
filter_order = 2;   % 4th order filter
cut_off = 4*(5/t_j(end)-t_j(1)); % cutoff frequency (Hz)
[num,den] = butter(filter_order,cut_off/(fs/2),'low');

v_j = cumtrapz(t_j,a_j); %m/s
d_j = cumtrapz(t_j,v_j); %m

figure;
set(gcf,'name','uncorrected velocity and displacement trajectories');
subplot(211)
hold on;
plot(t_j,v_j);
plot(t_j(ind_turn),v_j(ind_turn,:),'.k');
xlabel('time (s)'); ylabel('velocity (m/s)');
subplot(212)
hold on;
plot(t_j,cumtrapz(t_j,v_j),'color',[0.6,0.6,0.6]);
plot(t_j,d_j);
plot(t_j(ind_turn),d_j(ind_turn,:),'.k');
xlabel('time (s)'); ylabel('displacement (m)');


%% Define cone locations and shortest path through obstacle
c = zeros(7,3);
seg = 4.572; %15 ft in m
c(2,:) = [0,seg,0];
c(3,:) = [sqrt(seg^2 - (seg/2)^2),1.5*seg,0];
c(4,:) = [0,2*seg,0];
c(5,:) = [sqrt(seg^2 - (seg/2)^2),2.5*seg,0];
c(6,:) = [0,3*seg,0];
c(7,:) = [0,4*seg,0];

dc = interp1([t_j(1); t_j(ind_turn); t_j(end)],c,t_j,'linear');
   

%% Define and apply drift correction finding optimal heading correction
%Define initial, nominal heading angle
a0 = mean(a_j(1:round(0.5*ind_turn(1)),1:2));
a0 = a0./norm(a0);
theta0 = -atan2(a0(1),a0(2));

ind_still = find(ind);
f = @(x)vd_drift_correct(x,t_j,v_j,d_j,ind_turn,ind_still,c);

x0 = theta0;            
options = optimset('maxfunevals',10000,'maxiter',1000,'tolfun',1e-8,...
    'tolx',1e-8,'plotfcns',{@optimplotx,@optimplotresnorm});
lb = theta0-pi/2;
ub = theta0+pi/2;

theta = lsqnonlin(f,x0,lb,ub,options);

[vrc,drc] = apply_drift_correct(theta,t_j,v_j,d_j,ind_turn,ind_still,c);


%% Rotate data so initial angle aligned with cones
R = [cos(theta), sin(theta), 0;
    -sin(theta), cos(theta), 0;
        0,      0, 1];

dr = d_j * R.';
vr = v_j * R.';
ar = a_j * R.';
de = dr - dc;

arf = filtfilt(num,den,ar);


%% Plot drift corrected and aligned kinematics
figure;
set(gcf,'name','filtered accel in lab frame');
hold on;
plot(t_j,arf(:,1:2))
plot(t_j(ind_turn),arf(ind_turn,1:2),'ok','markerfacecolor','k');
xlabel('time (s)','fontsize',16); 
ylabel('acceleration (m/s^2)','fontsize',16);

figure;
set(gcf,'name','velocity drift error trajectories');
hold on;
plot(t_j(ind_turn),vr(ind_turn,:),'ok','markerfacecolor','k');
plot(t_j,vr(:,1:2)-vrc(:,1:2));
xlabel('time (s)','fontsize',16); 
ylabel('velocity drift (m/s)','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','corrected velocity trajectories');
hold on;
plot(t_j,vrc(:,1:2))
plot(t_j(ind_turn),vrc(ind_turn,1:2),'ok','markerfacecolor','k');
xlabel('time (s)','fontsize',16); 
ylabel('velocity (m/s)','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','corrected displacement trajectories');
hold on;
plot([t_j(ind_turn); t_j(end)],[dc(ind_turn,1:2); dc(end,1:2)],'ok','markerfacecolor','k');
plot(t_j,drc(:,1:2))
xlabel('time (s)','fontsize',16);
ylabel('displacement(m)','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','displacement drift error trajectories');
hold on;
plot([t_j(ind_turn); t_j(end)],[de(ind_turn,1:2); de(end,1:2)],'ok','markerfacecolor','k');
plot(t_j,dr(:,1:2)-drc(:,1:2));
xlabel('time (s)','fontsize',16); 
ylabel('displacement drift (m)','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','Corrected path');
hold on; grid on;
plot3(c(:,1),c(:,2),c(:,3),'ok');
plot3(dr(:,1),dr(:,2),dr(:,3),'m','linewidth',2);
plot3(drc(:,1),drc(:,2),drc(:,3),'r','linewidth',2);
plot3(drc(ind_turn,1),drc(ind_turn,2),drc(ind_turn,3),'ok','markerfacecolor','k');
plot3(dr(ind_turn,1),dr(ind_turn,2),dr(ind_turn,3),'ok','markerfacecolor','k');
xlabel('X (m)','fontsize',16); ylabel('Y (m)','fontsize',16);
set(gca,'fontsize',16);


%% Calculate additional kinamtic time series based on sacral trajectory
ap_dir = [ant_h, zeros(size(ant_h(:,1)))] * R.';
speed = sqrt(sum(vrc(:,1:2).^2,2)); %horizontal speed
speed_turn = speed(ind_turn); %horizontal

ind_mov = speed > 0.1;
tang = [vrc(:,1)./speed, ...
        vrc(:,2)./speed, ...
        zeros(size(vrc(:,3)))];

a_tang = dot(arf(ind_mov,:),tang(ind_mov,:),2);
a_resid = [arf(ind_mov,1:2), zeros(size(arf(ind_mov,3)))]-...
          [a_tang.*tang(ind_mov,1), a_tang.*tang(ind_mov,2), a_tang.*tang(ind_mov,3)];
norm_dir = [a_resid(:,1)./sqrt(sum(a_resid.^2,2)), ...
            a_resid(:,2)./sqrt(sum(a_resid.^2,2)),...
            a_resid(:,3)./sqrt(sum(a_resid.^2,2))];
a_norm = dot(arf(ind_mov,:),norm_dir,2);

ang_ap = ang_scan(atan2d(ap_dir(ind_mov,2),ap_dir(ind_mov,1)));
ang_tang = ang_scan(atan2d(tang(ind_mov,2),tang(ind_mov,1)));


[fit_pX, gof_pX] = fit(ang_tang,ang_ap,'poly1');
fit_ang = [fit_pX.p1, fit_pX.p2, gof_pX.rsquare, gof_pX.rmse]; %[slope,intercept, R^2, rms]


%% Plot more kinematics
figure;
set(gcf,'name','heading angle of ap and tang direction trajectories');
hold on;
plot(t_j(ind_mov),ang_ap,'b');
plot(t_j(ind_mov),ang_tang,'r');
plot(t_j(ind_turn),ang_ap(ind_turn-(find(ind_mov==1,1,'first')+1)),'ok','markerfacecolor','k');
xlabel('time (s)','fontsize',16); 
ylabel('heading angle (deg)','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','ap heading vs tang heading');
axis equal; grid on; hold on;
plot(ang_tang,ang_ap,'.b');
plot([min(ang_tang); max(ang_tang)],[min(ang_tang), 1; max(ang_tang), 1]*fit_ang(1:2).','k','linewidth',2);
plot([min(ang_tang); max(ang_tang)],[min(ang_tang), 1; max(ang_tang), 1]*fit_ang(1:2).'-fit_ang(4)*[1;1],'--k','linewidth',2);
plot([min(ang_tang); max(ang_tang)],[min(ang_tang), 1; max(ang_tang), 1]*fit_ang(1:2).'+fit_ang(4)*[1;1],'--k','linewidth',2);
xlabel('tangent heading (deg)','fontsize',16); 
ylabel('A-P heading (deg)','fontsize',16); 
set(gca,'fontsize',16);

figure;
set(gcf,'name','ap and tangential directions');
hold on;
plot(t_j(ind_mov),ap_dir(ind_mov,1:2))
plot(t_j(ind_mov),tang(ind_mov,1:2),'--')
plot(t_j(ind_turn),ap_dir(ind_turn,1:2),'ok','markerfacecolor','k');
xlabel('time (s)','fontsize',16); 
ylabel('direction','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','normal and tangential accelerations');
hold on; grid on;
plot(t_j(ind_mov),a_tang)
plot(t_j(ind_mov),a_norm,'r')
legend('tangential','normal');
plot(t_j(ind_turn),a_tang(ind_turn-(find(ind_mov==1,1,'first')+1)),'ok','markerfacecolor','k');
plot(t_j(ind_turn),a_norm(ind_turn-(find(ind_mov==1,1,'first')+1)),'ok','markerfacecolor','k');
xlabel('Time (s)','fontsize',16); 
ylabel('Tangential Acceleration (m/s^2)','fontsize',16);
set(gca,'fontsize',16);

figure;
set(gcf,'name','horizontal speed');
hold on;
plot(t_j,speed)
plot(t_j(ind_turn),speed(ind_turn),'ok','markerfacecolor','k');
xlabel('time (s)','fontsize',16); 
ylabel('horizontal speed (m/s)','fontsize',16);
set(gca,'fontsize',16);