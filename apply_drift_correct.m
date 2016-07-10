function [vc,dc] = apply_drift_correct(x,t_j,v,d,ind,ind_still,c)

%Apply rotation correction
R = [cos(x), sin(x), 0;
    -sin(x), cos(x), 0;
          0,      0, 1];

vr = v * R.';
dr = d * R.';

%Define and apply correction
t=t_j-t_j(1);
wi = [1; 1; 1; 1; 1];
ti = t(ind);
di = dr(ind,:)-c(2:6,:);

%Initial zero velocity constraint, adjusted t0, cubic
i0 = ind_still(find(ind_still<ind(1),1,'last'));
t0 = t(i0); d0 = dr(i0,:); v0 = vr(i0,:);
t6 = t(end); d6 = dr(end,:);

T = [1, t0, t0^2, t0^3, 0, 0, 0;
    1, t6, t6^2, t6^3, 0, 0, 0;
    0, 1, 2*t0, 3*t0^2, 0, 0, 0;
    2*wi.'*ti.^0, 2*wi.'*ti.^1, 2*wi.'*ti.^2, 2*wi.'*ti.^3, -1, -1, 0;
    2*wi.'*ti.^1, 2*wi.'*ti.^2, 2*wi.'*ti.^3, 2*wi.'*ti.^4, -t0, -t6, -1;
    2*wi.'*ti.^2, 2*wi.'*ti.^3, 2*wi.'*ti.^4, 2*wi.'*ti.^5, -t0^2, -t6^2, -2*t0;
    2*wi.'*ti.^3, 2*wi.'*ti.^4, 2*wi.'*ti.^5, 2*wi.'*ti.^6, -t0^3, -t6^3, -3*t0^2];

y(:,1) = [d0(1);
    d6(1)-c(end,1);
    v0(1);
    2*wi.'*di(:,1);
    2*(wi.*ti).'*di(:,1);
    2*(wi.*ti.^2).'*di(:,1);
    2*(wi.*ti.^3).'*di(:,1)];
y(:,2) = [d0(2);
    d6(2)-c(end,2);
    v0(2);
    2*wi.'*di(:,2);
    2*(wi.*ti).'*di(:,2);
    2*(wi.*ti.^2).'*di(:,2);
    2*(wi.*ti.^3).'*di(:,2)];
y(:,3) = [d0(3);
    d6(3)-c(end,3);
    v0(3);
    2*wi.'*di(:,3);
    2*(wi.*ti).'*di(:,3);
    2*(wi.*ti.^2).'*di(:,3);
    2*(wi.*ti.^3).'*di(:,3)];

coe = T \ y;

%Cubic w/ i0
ind_s = ind_still(ind_still<=i0);
vc(1:i0,:) = vr(1:i0,:) - interp1(t(ind_s),vr(ind_s,:),t(1:i0),'linear');
dc(1:i0,:) = cumtrapz(t(1:i0),vr(1:i0,:));

t = t(i0:end);
len = length(dr);
dc(i0:len,:) = dr(i0:end,:) - [t.^0, t, t.^2, t.^3]*coe(1:4,:) + (t.^0)*dc(i0,:);
            
vc(i0:len,:) = vr(i0:end,:) - [t.^0, 2*t, 3*t.^2]*coe(2:4,:);

end

