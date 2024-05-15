function [r, theta] = getPos(tdtx1, tdtx2,dist)
% [r, theta] = getPos(tdtx1, tdtx2,dist)
% Finds the position of the object in x-y space based of
% the delays for each RADAR transceiver.

    c = 299792458;
    c1 = tdtx1*c;
    c2 = tdtx2*c;
    A = [2*c1 -1*dist; 2*c2 dist];
    b = [(c1^2-dist/4); (c2^2-dist/4)];
    ry = A\b;
    r = ry(1);
    theta = asin(ry(2)/ry(1));
end