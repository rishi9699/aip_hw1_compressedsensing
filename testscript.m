num_rows = 120; 
num_columns = 210; % Using a 120x210 window
num_frames=3;

video=mmread('./HW1/cars.avi');
frames = zeros(num_rows, num_columns, num_frames);
for i=1:num_frames
    frames(:,:,i) = rgb2gray(video.frames(i).cdata(169:288, 143:352,:));
end

random_pattern = binornd(1,0.7,num_rows,num_columns,num_frames);

% Computing coded snapshot
coded_snapshot = sum(frames.*random_pattern, 3) %+ 2*randn(num_rows,num_columns); % Adding Gaussian Random Noise
imshow(uint8(coded_snapshot/num_frames)); % Display coded snapshot


patch_size=8;
block_size = patch_size^2;

%% Constructing the 3D DCT psi matrix
psi = zeros(block_size*num_frames, block_size*num_frames);

const_coeff = 1/sqrt(num_frames*patch_size*patch_size);
coeff =ones(1, 3);

for w=0:num_frames-1
    if w==0
        coeff(1) = 1;
    else
        coeff(1) = sqrt(2);
    end
    for v=0:patch_size-1
        if v==0
            coeff(2) = 1;
        else
            coeff(2) = sqrt(2);
        end
        for u=0:patch_size-1
            if u==0
                coeff(3) = 1;
            else
                coeff(3) = sqrt(2);
            end
            
            row_number = block_size*w + patch_size*v + (u+1);
            
            
            for t=0:num_frames-1
                for m=0:patch_size-1
                    temp = block_size*t + patch_size*m;
                    for n=0:patch_size-1
                        psi(row_number, temp + n + 1) = const_coeff * prod(coeff) * cos((pi*(2*n+1)*u)/(2*patch_size)) * cos((pi*(2*m+1)*v)/(2*patch_size)) * cos((pi*(2*t+1)*w)/(2*num_frames));
                    end
                end
            end
            
            
        end
    end
end

%%
reconstructed_frames = zeros(num_rows, num_columns, num_frames);
phi = zeros(block_size, block_size*num_frames);
all_weights = zeros(num_rows, num_columns);

for patch_start_row = 1:(num_rows-patch_size+1) % Iterating over all possible patches
    for patch_start_column = 1:(num_columns-patch_size+1)
        

        % Constructing the phi matrix
        for i=0:num_frames-1
            phi(:, (i*block_size+1):(i+1)*block_size) = diag(reshape(random_pattern(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1,i+1).',[1, block_size]));
        end

        % Implementing OMP
        indices = zeros(1, block_size);
        support_set = zeros(block_size, block_size);
        
        % Check OMP algorithm ?
        A=phi*psi;
        normalized_A = normc(A);
        measurement = reshape(coded_snapshot(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1).', [1,block_size]).';
        residual = measurement;
        
        iter=1;
        norm_residual=2305;
        while norm_residual>600 && iter<(block_size+1)
        %for iter=1:block_size % Running the OMP algorithm for a maximum of block_size iterations for well-defined pseudo inverse
            [~, i] = max(abs(residual.'*normalized_A));
            indices(iter) = i;
            support_set(:,iter) = A(:, i);
            theta = pinv( support_set(:,1:iter) ) * measurement;
            residual = measurement - support_set(:,1:iter) * theta;
            norm_residual = norm(residual);
            iter = iter + 1;
        %end
        end
        
        iter = iter - 1;
        % Converting signal from DCT basis to Dirac basis
        x = zeros(block_size*num_frames, 1);
        for pos=1:iter
            x(indices(pos)) = theta(pos);
        end
        %f = psi*x; % DCT forward formula ?
        f = sum(psi.*x);
        f=reshape(f, [patch_size, patch_size, num_frames]); % Patch frames in the Dirac basis
        f = permute(f, [2 1 3]);

        % Rejoining the patches
        weights = all_weights(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1);

        reconstructed_frames(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1,:) ...
         = (reconstructed_frames(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1,:) .* (weights./(weights+1))) + (f ./ (weights+1));

        all_weights(patch_start_row:patch_start_row+patch_size-1,patch_start_column:patch_start_column+patch_size-1) ...
            = (weights+1);
    end
    disp(patch_start_row)
end

for p=1:num_frames % Displaying the recostructed frames
    figure
    imshow(uint8(reconstructed_frames(:,:,p)))
end

% % Display RRMSE
% for p=1:num_frames
%     sqrt(sum((uint8(reconstructed_frames(:,:,p)) - uint8(frames(:,:,p))).^2))/sqrt(sum(uint8(frames(:,:,p)).^2))
% end
