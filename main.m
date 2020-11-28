%%% Predicts part quality considering locating surfaces and machining error%%%

clear; clc; close all;
include_namespace_dq

% The pos of 5 locating surfaces and tool deviations
locatingSurfs = getLocatingSurfs();
f = [1:3; 4:6; 7:9; 10:12; 13:15]';
toolDeviations = [0, -0.01, 0.02, -0.03, 0.04];

ParallelismS1s = [];
ParallelismS8s = [];
for e = 1:5 % 5 experiments
    disp(strcat('Experiment ', num2str(e)))

    % get fixture and workpiece info
    fixDatums = locatingSurfs(:, f(:, e)');
    v8WorkPiece = getWorkpiece2();

    % Run operations
    %%
    R = zeros(8, 1); R(1) = 1;
    v8WorkPiece = operationStn1(v8WorkPiece, fixDatums, R, toolDeviations(e));
    v8WorkPiece = operationStn2(v8WorkPiece, fixDatums, R, toolDeviations(e));

    % final check primary
    %%
    tol = 10^-6;
    err = getProjectedTestPointsPlucker2(fixDatums(:, 1), [0, 0, 0], fixDatums(2:4, 1)) - ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 3), [0, 0, 0], fixDatums(2:4, 1));

    % Move to insection station and conduct measurement
    %%
    v8WorkPiece = moveToInspectionStation(v8WorkPiece, 3, 7);

    S1 = [getProjectedTestPointsPlucker2(v8WorkPiece(:, 2), [0, 0, 120], [0, 0, 1]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 2), [0, 200, 120], [0, 0, 1]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 2), [200, 0, 120], [0, 0, 1]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 2), [200, 200, 120], [0, 0, 1])];

    ParallelismS1s = [ParallelismS1s; round(max(S1(:, 3))-min(S1(:, 3)), 3)];

    S8_1 = [getProjectedTestPointsPlucker2(v8WorkPiece(:, 9), [40, 160, 100], [0, 1, 0]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 9), [40, 160, 120], [0, 1, 0]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 9), [160, 160, 100], [0, 1, 0]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 9), [160, 160, 120], [0, 1, 0])];

    S8_2 = [getProjectedTestPointsPlucker2(v8WorkPiece(:, 5), [40, 160, 100], [0, 1, 0]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 5), [40, 160, 120], [0, 1, 0]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 5), [160, 160, 100], [0, 1, 0]); ...
        getProjectedTestPointsPlucker2(v8WorkPiece(:, 5), [160, 160, 120], [0, 1, 0])];

    S8_1_2 = abs(S8_1-S8_2);
    ParallelismS8s = [ParallelismS8s; round(max(S8_1_2(:, 2))-min(S8_1_2(:, 2)), 3)];
end
%
% All predicted parallelism values
Parallelisms = [ParallelismS1s, ParallelismS8s]
%
%
%%
function v8WorkPiece = operationStn1(v8WorkPiece, fixDatums, RMach, ToolDeviation)
disp('Running Station 1 ...')

include_namespace_dq

% assemble to fixture
stn = 1;
v8WorkPiece = assemblePrimSecTert(v8WorkPiece, fixDatums, stn);

% machining error
s1CuttingPlane = vec8(k_+E_*120);
TMach = vec8(1-E_*0.5*ToolDeviation*k_);
RT = DQmult(TMach, RMach);

% replace with the machined surf
v8WorkPiece(:, 3) = moveFeaturesBy(s1CuttingPlane, RT);

end
%%
function v8WorkPiece = operationStn2(v8WorkPiece, fixDatums, RMach, ToolDeviation)
disp('Running Station 2 ...')

include_namespace_dq

% rotate first
v8WorkPiece = RotatePart180DegXaxis(fixDatums, v8WorkPiece, 2);
v8WorkPiece = flipNormalsWorkpiece(v8WorkPiece);

% Assemble to fixture
stn = 2;
v8WorkPiece = assemblePrimSecTert(v8WorkPiece, fixDatums, stn);

% Nominal cutting planes
s2CuttingPlaneH = vec8(k_+E_*100);
s2CuttingPlaneV1 = vec8(j_+E_*160);

% machining error
TMach = vec8(1-E_*0.5*ToolDeviation*k_);
TMachV = vec8(1-E_*0.5*ToolDeviation*j_);
RT = DQmult(TMach, RMach);
RTV = DQmult(TMachV, RMach);

v8WorkPiece(:, 2) = moveFeaturesBy(s2CuttingPlaneH, RT);
v8WorkPiece(:, 9) = moveFeaturesBy(s2CuttingPlaneV1, RTV);
end
%%
%
function v8WorkPiece = assemblePrimSecTert(v8WorkPiece, fixDatums, stn)
% Assembles to primary, secondary, and tertiary locating surfaces
include_namespace_dq
tol = 10^-6;

% get intersection lines and Direction based on fixture orientation
ULPS = checkAndFlipNormalDir(normalize(cross(fixDatums(2:4, 1), fixDatums(2:4, 2))));
ULPT = checkAndFlipNormalDir(normalize(cross(fixDatums(2:4, 1), fixDatums(2:4, 3))));
assert(ULPS(1) > 0)
assert(ULPT(2) > 0)

% Datum
switch stn
    case 1
        pDatum = 1;
        sDatum = 5;
    case 2
        pDatum = 3;
        sDatum = 7;
end

% Assembly of primary feature
%%
[R, errD] = getErrTransformationParas(fixDatums(:, 1), v8WorkPiece(:, pDatum), [], [], 3);
T = vec8(1-E_*errD/2*fixDatums(2:4, 1));
RT = DQmult(T, R);

v8WorkPiece = moveFeaturesBy(v8WorkPiece, RT);
err = getProjectedTestPointsPlucker2(fixDatums(:, 1), [0, 0, 0], fixDatums(2:4, 1)) - ...
    getProjectedTestPointsPlucker2(v8WorkPiece(:, pDatum), [0, 0, 0], fixDatums(2:4, 1));
assert(norm(err) < tol);
disp('Passed primary check')

% Assembly of secondary feature
%%
[pts45T, pts45B] = getPointsOnSecTerDatumShortFix(v8WorkPiece, fixDatums, sDatum, 's'); %Project on Fixture and rotate
[p4Bfix, p5Bfix, p4Tfix, p5Tfix] = projectOntoSecondaryFix(pts45B, pts45T, fixDatums, ULPT);

% Rotate part first to make it parallel with secondary locator
p4Tp5Tvec = pts45T(1, :) - pts45T(2, :); % vector connecting both points
p4Tp5TvecFix = p4Tfix' - p5Tfix';

[R, errD] = getErrTransformationParas(normalize(p4Tp5TvecFix), normalize(p4Tp5Tvec), [], [], 3);
v8WorkPiece = moveFeaturesBy(v8WorkPiece, R);

% check primary again
err = getProjectedTestPointsPlucker2(fixDatums(:, 1), [0, 0, 0], fixDatums(2:4, 1)) - ...
    getProjectedTestPointsPlucker2(v8WorkPiece(:, pDatum), [0, 0, 0], fixDatums(2:4, 1));
assert(norm(err) < tol);

% check for equal distance from locators after rotation
%%
[pts45T, pts45B] = getPointsOnSecTerDatumShortFix(v8WorkPiece, fixDatums, sDatum, 's');
[p4Bfix, p5Bfix, p4Tfix, p5Tfix, t4B, t5B, t4T, t5T] = ...
    projectOntoSecondaryFix(pts45B, pts45T, fixDatums, ULPT);
errT = [p4Tfix'; p5Tfix'] - pts45T(1:2, :);
errB = [p4Bfix'; p5Bfix'] - pts45B(1:2, :);

% assert(errT(1, 2) - errT(2, 2) < tol); % checks equal distance

% take the shortest of the two pairs
if abs(errT(1, 2)) < abs(errB(2, 2))
    err = errT(1, :);
    disp('Contact at the top')
else
    err = errB(1, :);
    disp('Contact at the bottom')
end
T = vec8(1-E_*err/2);
v8WorkPiece = moveFeaturesBy(v8WorkPiece, T);

% check primary
err = getProjectedTestPointsPlucker2(fixDatums(:, 1), [0, 0, 0], fixDatums(2:4, 1)) - ...
    getProjectedTestPointsPlucker2(v8WorkPiece(:, pDatum), [0, 0, 0], fixDatums(2:4, 1));
assert(norm(err) < tol);

% check secondary
[pts45T, pts45B] = getPointsOnSecTerDatumShortFix(v8WorkPiece, fixDatums, sDatum, 's');
errT = [p4Tfix'; p5Tfix'] - pts45T;
errB = [p4Bfix'; p5Bfix'] - pts45B;

try
    errt = errT(1, :);
    assert(norm(errT(1, :)) < tol)
catch
    errb = errB(1, :);
    assert(norm(errB(1, :)) < tol)
end
disp('Passed secondary check')

% Assembly of tertiary feature
%%
tDatum = 6;
[pts45Tt, pts45Bt] = getPointsOnSecTerDatumShortFix(v8WorkPiece, fixDatums, tDatum, 't');
vTFeature = [pts45Tt; pts45Bt];

R = zeros(8, 1); R(1) = 1;
err = projOnToFixtureGetErr(fixDatums(2:5, 3), vTFeature, ULPS);
T = vec8(1-E_*err*ULPS/2);
v8WorkPiece = moveFeaturesBy(v8WorkPiece, T);

% check contact with primary again
err = getProjectedTestPointsPlucker2(fixDatums(:, 1), [0, 0, 0], fixDatums(2:4, 1)) - ...
    getProjectedTestPointsPlucker2(v8WorkPiece(:, pDatum), [0, 0, 0], fixDatums(2:4, 1));
assert(norm(err) < tol);
%
disp('Passed tertiary check')
end
%%
function err = projOnToFixtureGetErr(fixPlane, vFeature, dir)
iPts = [];
dists = [];

for v = 1:size(vFeature, 1)
    L = Plucker.pointdir(vFeature(v, :), dir);
    [p, t] = L.intersect_plane2(fixPlane); % intersection to Fix
    if ~isempty(p)
        dist = norm([p' - vFeature(v, :)]);
        iPts = [iPts; dist, p'];
    end
end

% get the shortest distance from four vertices
err = min(abs(iPts(:, 1)));
shortestDist = iPts(find(min(abs(iPts(:, 1))) == abs(iPts(:, 1))), 2:4);
shortestDistNorm = norm(shortestDist);
end
%%
function [pts45T, pts45B] = getPointsOnSecTerDatumShortFix(v8WorkPiece, fixDatums, vDatum, dtm)
include_namespace_dq

% two vertical lines (% vDatum vertical datum)
if dtm == 's'
    LS5 = Plucker.planes(v8WorkPiece(2:5, 6), v8WorkPiece(2:5, vDatum));
    LS7 = Plucker.planes(v8WorkPiece(2:5, 8), v8WorkPiece(2:5, vDatum));

else dtm == 't';
    LS5 = Plucker.planes(v8WorkPiece(2:5, 5), v8WorkPiece(2:5, vDatum));
    LS7 = Plucker.planes(v8WorkPiece(2:5, 7), v8WorkPiece(2:5, vDatum));
end
% Get a plane parallel to the primary passing through height 75 mm
plane75 = vec8(getDQPlane(fixDatums(2:4, 1), [0, 0, 75]));

p4T = LS5.intersect_plane3(plane75(2:5))';
p5T = LS7.intersect_plane3(plane75(2:5))';
p4B = LS5.intersect_plane3(fixDatums(2:5, 1))';
p5B = LS7.intersect_plane3(fixDatums(2:5, 1))';

pts45T = [p4T; p5T];
pts45B = [p4B; p5B];
end
%%
function[p4Bfix, p5Bfix, p4Tfix, p5Tfix, t4B, t5B, t4T, t5T] = projectOntoSecondaryFix(pts45B, pts45T, fixDatums, ULPT)
L4B = Plucker.pointdir(pts45B(1, :), ULPT);
L5B = Plucker.pointdir(pts45B(2, :), ULPT);
L4T = Plucker.pointdir(pts45T(1, :), ULPT);
L5T = Plucker.pointdir(pts45T(2, :), ULPT);

[p4Bfix, t4B] = L4B.intersect_plane2(fixDatums(2:5, 2));
[p5Bfix, t5B] = L5B.intersect_plane2(fixDatums(2:5, 2));
[p4Tfix, t4T] = L4T.intersect_plane2(fixDatums(2:5, 2));
[p5Tfix, t5T] = L5T.intersect_plane2(fixDatums(2:5, 2));
end
