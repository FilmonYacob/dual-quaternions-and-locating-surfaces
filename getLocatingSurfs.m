function locatingSurfs = getLocatingSurfs()
% Inputs of experiments 1-5
axis = [0, 1, 0; 0, 0, 1; 0, 0, 1];
thetas = [[0.3; -0.1; 0.1], [-0.1; 0.4; 0.2], [-0.2; -0.2; 0.3], [0.2; 0.5; 0.2], [-0.03; -0.02; 0.2]];
ds = [[0, 0, 0.1, 0, 0.05, 0, 0.5, 0, 0]; [0, 0, 0.4, 0, 0.05, 0, 0.1, 0, 0]; [0, 0, 0.2, 0, 0.2, 0, -0.3, 0, 0]; [0, 0, 0.1, 0, 0.3, 0, -0.2, 0, 0]; [0, 0, 0.1, 0, 0.1, 0, 0.1, 0, 0]]';

locatingSurfs = [];

include_namespace_dq
for cols = 1:5
    theta = thetas(:, cols);
    d = [ds(1:3, cols)'; ds(4:6, cols)'; ds(7:9, cols)'];
    for i = 1:3
        R = normalize(rot2dquat(theta(i), axis(i, :)));
        T = vec8(1-E_*d(i, 1:3)/2);
        RT = DQmult(T, R);
        if i == 1
            plane_ = vec8(k_);
        elseif i == 2
            plane_ = vec8(j_);
        else
            plane_ = vec8(i_);
        end
        
        k = DQ(moveFeaturesBy(plane_, RT))
        locatingSurfs = [locatingSurfs, moveFeaturesBy(plane_, RT)];
    end
end
end