function e = vd_drift_correct(x,t,v,d,ind,ind_still,c)
[~,dc] = apply_drift_correct(x,t,v,d,ind,ind_still,c);

%Minimize distance from cones
e = [c(2:6,1) - dc(ind,1); c(2:6,2) - dc(ind,2)];
end