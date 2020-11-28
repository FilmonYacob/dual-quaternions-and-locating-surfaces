function v8WorkPiece = RotatePart180DegXaxis(fixDatums, v8WorkPiece, stn)
include_namespace_dq

R = rot2dquat(180, [1, 0, 0]');
v8WorkPiece = moveFeaturesBy(v8WorkPiece, R);
trans1 = 1 - E_ * j_ / 2;
trans2 = 1 - E_ * 45 * k_ / 2;
T = DQmult(trans1.vec8,trans2.vec8);
v8WorkPiece = moveFeaturesBy(v8WorkPiece, T);

% check for large deviations
 switch stn
    case 1
        sDatum = 5;
    case 2
        sDatum = 7;
end
sd = getProjectedTestPointsPlucker2(fixDatums(:,2), [0, 0, 0], fixDatums(2:4,2)) ...
    -getProjectedTestPointsPlucker2(v8WorkPiece(:, sDatum), [0, 0, 0], fixDatums(2:4,2));
if v8WorkPiece(5,sDatum)*v8WorkPiece(3,sDatum) < 0
    T = 1 - E_ * sd / 2;
    v8WorkPiece = moveFeaturesBy(v8WorkPiece, T.vec8); 
end
end