function xnew = ang_scan(x)
%Function to fix jumps in angle data, expects angle data in degrees

%Calc change in angle from sample to sample
dx = diff(x);

%Identify if there are jumps in angle data greater than 180 deg
jps = find(abs(dx)>180);

if isempty(jps) %No jumps
    xnew = x;
else %Jumps
    xnew = x;
    jp_vals = dx(jps);
    i=1;
    jump_cum = 0;
    while i<=length(jp_vals)
        jump_start = jps(i)+1;
        jump_dir = sign(jp_vals(i));
        try
            jump_end = jps(i+1);
            xnew(jump_start:jump_end) = x(jump_start:jump_end)-jump_dir*360+jump_cum;
            i = i + 1;
            jump_cum = jump_cum-jump_dir*360;
        catch
            jump_end = length(x);
            xnew(jump_start:jump_end) = x(jump_start:jump_end)-jump_dir*360+jump_cum;
            i = length(jp_vals)+1;
            jump_cum = jump_cum-jump_dir*360;
            xnew(jump_end:end) = x(jump_end:end)+jump_cum;
        end
    end
end
end