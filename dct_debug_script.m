testm = ones(1,192);
u=0;
v=0;
w=0;
for t=0:num_frames-1
    for m=0:patch_size-1
        temp = block_size*(t) + patch_size*(m);
        for n=0:patch_size-1
            disp(t+" "+m+" "+n)
            disp(cos((pi*(2*n+1)*u)/(2*patch_size)));
            disp(cos((pi*(2*m+1)*v)/(2*patch_size)));
            disp(cos((pi*(2*t+1)*w)/(2*num_frames)));
            testm(temp + n + 1) = sqrt(1/192) * cos((pi*(2*n+1)*u)/(2*patch_size)) * cos((pi*(2*m+1)*v)/(2*patch_size)) * cos((pi*(2*t+1)*w)/(2*num_frames));
        end
    end
end

for t=0:num_frames-1
    for m=0:patch_size-1
        for n=0:patch_size-1
            disp(t+" "+m+" "+n)
            disp(cos((pi*(2*n+1)*u)/(2*patch_size)));
            disp(cos((pi*(2*m+1)*v)/(2*patch_size)));
            disp(cos((pi*(2*t+1)*w)/(2*num_frames)));
            testm((t+1)*(m+1)*(n+1)) = sqrt(1/192) * cos((pi*(2*n+1)*u)/(2*patch_size)) * cos((pi*(2*m+1)*v)/(2*patch_size)) * cos((pi*(2*t+1)*w)/(2*num_frames));
        end
    end
end
