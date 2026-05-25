function trajNorm = normalizeTrajectory(traj)
% Normalize trajectory by path length for comparison

pathLength = sum(sqrt(sum(diff(traj, 1, 1).^2, 2)));
if pathLength > 0
    trajNorm = traj / pathLength;
else
    trajNorm = traj;
end
end