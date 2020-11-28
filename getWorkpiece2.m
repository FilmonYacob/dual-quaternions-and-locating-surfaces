function v8WorkPiece = getWorkpiece2()
% creates all features of the part

include_namespace_dq

S0 = k_;
S1 = k_ + E_ * 100;
S2 = k_ + E_ * 120;
S3 = i_ + E_ * 40;
S4 = j_ + E_ * 0;
S5 = i_ + E_ * 0;
S6 = j_ + E_ * 200;
S7 = i_ + E_ * 200;
S8 = j_ + E_ * 160;

v8WorkPiece = [S0.vec8, S1.vec8, S2.vec8, S3.vec8, S4.vec8, ...
    S5.vec8, S6.vec8, S7.vec8, S8.vec8];
end