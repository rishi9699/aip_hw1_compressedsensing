num_rows=288;
num_columns=352;
num_frames=3;

video=mmread('../HW1/cars.avi');
frames = uint8(zeros(num_rows, num_columns, num_frames));
for i=1:num_frames
    frames(:,:,i) = im2gray(video.frames(i).cdata);
end

figure
imshow(frames(:,:,1))
figure
imshow(frames(:,:,2))
figure
imshow(frames(:,:,3))

random_pattern = uint8(binornd(1,0.5,num_rows,num_columns,num_frames));

% Computing coded snapshot
summand_matrices = frames.*random_pattern;
coded_snapshot = uint8(sum(summand_matrices, 3));
imshow(coded_snapshot);


%%

patch_size=8;
block_size = patch_size^2;

% Constructing the 3D DCT psi matrix
psi = zeros(block_size*num_frames, block_size*num_frames);

coeff=1/sqrt(num_frames*path_size*patch_size);
for w=0:num_frames-1
    if w==0
        %nothing
    else
        coeff = coeff * sqrt(2);
    end
    for v=0:patch_size-1
        if v==0
            %nothing
        else
            coeff = coeff * sqrt(2);
        end
        for u=0:patch_size-1
            if u==0
                %nothing
            else
                coeff = coeff * sqrt(2);
            end
            row_number = (w+1 * v+1 * u+1);
            
            for t=0:num_frames-1
                for m=0:patch_size-1
                    temp = (t+1 * m+1);
                    for n=0:patch_size-1
                        psi(row_number, temp*(n+1)) = coeff * cos((pi*(2*n+1)*u)/patch_size) * cos((pi*(2*m+1)*v)/patch_size) * cos((pi*(2*t+1)*w)/num_frames);
                    end
                end
            end
        end
    end
end
            coeff = 1/sqrt(num_frames*path_size*patch_size);


patch_start_row = 10;
patch_start_column = 10;

% Constructing the phi matrix
phi = zeros(block_size, block_size*num_frames);
for i=0:num_frames-1
    phi(:, (i*block_size+1):(i+1)*block_size) = diag(reshape(random_pattern(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1,i+1),[block_size,1]));
end
