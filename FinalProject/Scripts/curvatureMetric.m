% Behavioral: Find a curvature metric (ratio of actual path length to 
% straight line distance between start and end). Plot the distribution of 
% curvature with a bootstrapped confidence interval for curved and straight 
% path conditions. This would be a behavioral verification that the animal 
% actually did the curved paths and not just the straight length ones. 

num_trials = length(R);
range = [-200, 800];
num_directions = 36;
mazeCategory = [R.mazeCategory];
moveOnsetTimes = [R.moveOnsetTime];
straight_idx = mazeCategory == 'noBarriers';
curved_idx   = mazeCategory == 'barriers';

trial = 1;
hand_times = R(trial).HAND.times;
x_pos = R(trial).HAND.X;
y_pos = R(trial).HAND.Y;
moveOnsetTime = moveOnsetTimes(trial);
startTime = moveOnsetTime + range(1);
endTime = moveOnsetTime + range(2);

[~, startTime_idx] = min(abs(hand_times - startTime)); 
[~, moveOnsetTime_idx] = min(abs(hand_times - moveOnsetTime));
[~, endTime_idx] = min(abs(hand_times - endTime));

figure;

