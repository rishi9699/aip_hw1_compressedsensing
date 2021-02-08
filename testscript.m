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

%Computing coded snapshot
summand_matrices = frames.*random_pattern;
coded_snapshot = uint8(sum(summand_matrices, 3));
imshow(coded_snapshot);


%%
patch_size=8;
block_size = patch_size^2;

patch_start_row = 10;
patch_start_column = 10;

%Constructing the A matrix
A = zeros(block_size, block_size*num_frames);
for i=0:num_frames-1
    A(:, (i*block_size+1):(i+1)*block_size) = diag(reshape(random_pattern(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1,i+1),[block_size,1]));
end

