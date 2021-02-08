A = rand(4, 12);
normalized_A = normc(A);
measurement = randn(1, 4).';
block_size = 4;
indices = zeros(1, block_size);
support_set = zeros(block_size, block_size);

residual = measurement;
for iter=1:block_size
    [~, i] = max(abs(residual.'*normalized_A));
    indices(iter) = i;
    support_set(:,iter) = A(:, i);
    theta = pinv( support_set(:,1:iter) ) * measurement;
    residual = measurement - support_set(:,1:iter) * theta;
end

x = zeros(block_size*num_frames, 1);
for pos=1:block_size
    x(indices(pos)) = theta(pos);
end
    
%Although the sparse theta vector obtained on multiplying with A gives ~y,
%how do we know it is the actual solution?
