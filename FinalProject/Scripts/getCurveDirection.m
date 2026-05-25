function curveDir = getCurveDirection(R)
% Classify each trial as clockwise, counterclockwise, or straight.
% Straight vs curved comes from the maze category. For curved trials, the
% sign of the path's deviation from the straight line determines CW vs CCW.


nTrials = length(R);
curveDir = strings(nTrials, 1);
isStraight = [R.mazeCategory] == 'noBarriers';

for tr = 1:nTrials
    if isStraight(tr)
        curveDir(tr) = "Straight";
        continue
    end

    % Curved trial
    hx = R(tr).HAND.X;
    hy = R(tr).HAND.Y;
    ht = R(tr).HAND.times - R(tr).moveOnsetTime;
    movDur = R(tr).moveEndsTime - R(tr).moveOnsetTime;
    mask = ht >= 0 & ht <= movDur;

    x = hx(mask) - hx(find(mask, 1));
    y = hy(mask) - hy(find(mask, 1));

    % Straight line
    targetDir = [x(end), y(end)] / norm([x(end), y(end)]);

    % Trajectory relative to straight line
    midPt = [x(round(end/2)), y(round(end/2))];
    sideValue = targetDir(1) * midPt(2) - targetDir(2) * midPt(1);

    if sideValue > 0
        curveDir(tr) = "CCW";
    else
        curveDir(tr) = "CW";
    end
end
end