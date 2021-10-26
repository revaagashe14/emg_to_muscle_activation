%leastSquares(0,0,0,0,0,0,0,0,0,0,0)
fullexample
function f = leastSquares(input)
    
    A1 = input(1);    
    A2 = input(2);
    A3 = input(3);
    A4 = input(4);
    A5 = input(5);
    A6 = input(6);
    A7 = input(7);
    A8 = input(8);
    din = input(9);
    y1in = input(10);
	y2in = input(11);
    Ain = [A1; A2; A3; A4; A5; A6; A7; A8];
    
    t_end = 4000;
% 
    b1 = y1in + y2in;
    b2 = y1in * y2in;
    alpha1 = 1 + b1 + b2;
  
    f = 0;
    %din = din/5;
    dataMatrix = load('Data/s1.mat', 'dsfilt_emg', 'finger_kinematics'); 
    %dataMatrix
    for trial = 1:1 %should be 1:3
        for activity = 5:5 %should be 1:7
            
            emg_signal = dataMatrix.('dsfilt_emg');
            trial_activity_data = emg_signal{trial, activity};

            kinematics = dataMatrix.('finger_kinematics');
            kinematics_data = kinematics{trial, activity};
            
            u = zeros(4000, 8); % neural activation (used to calc v)
            v = zeros(4000, 8); % muscle activation
            N = 15; %should be 15 for 15 dofs for the orig dataset
            theta_targets = zeros(4000, 15);
            theta_estimates = zeros(4000, 15);  %use gpr to find this, compute/estimate from EMG signal data
  
            for t = 1:4000  %t is the timestep (out of 4000)              
                    %calculate theta_targets = correct joint angles for a
                    %given joint positions for a specific trial and activity and timestep
                    theta_targets(t,:) = calculateExpectedJointAngles(t, kinematics_data);
                    % theta_targets
                    
                   % trial_activity_data
                   e = zeros(8, 4000);
                % calculating e below
                    for x = 1:8  %x is the muscle (out of 8)
                        %rectify using absolute value
                        for i = 1:4000
                            e(x,i) = abs(trial_activity_data(i, x));
                        end
                        % fix(trial_activity_data(:, x));  % change this to get it from data
                        % normalize
                        e(x, :) = e(x, :)/max(e(x, :));
                                              
                        % filter using a 2nd order butterworth filter
                        fc = 4;
                        fs = 200;

                        [b,a] = butter(2,fc/(fs/2));

                        e(x, :) = filter(b,a,e(x, :));
                        %evaluate = e(t)
                        
                        %make e continuous so t-din can be evaluated
                        % funcEx = ts2func(e(x,:));
                        
                        % find u(t, x): (later account for base cases: what if din>t? are the if cases right)
                        t
                        din
                        alpha1;
                        % fval = funcEx(t-din);

                        if (t*5-din >=3) && (t*5-din <= 20000)
                            % u(t, x) = alpha1*funcEx(t-din);
                            u(t,x) = alpha1*e(x, round((t*5-din)/5));
                            fprintf("here");
                            
                        end
                        if (t-1 >=1)
                            u(t, x) = u(t,x) - b1*u(t-1, x);
                            fprintf("here2");
                        end

                        if (t-2 >=1)
                            u(t, x) = u(t,x) - b2 * u(t-2, x);
                            fprintf("here3");
                        end  
                        test = u(t,x);

                        % calculate v(t,x)
                        v(t,x) = (exp(Ain(x)* u(t,x)) - 1) / (exp(Ain(x))-1);    
                        test2 = v(t,x);
                        
                    end
                    %center v(t,:) to have zero mean and unit variance
                    %v(t,:) = (v(t,:) - mean(v(t,:)))/std(v(t,:));
                                      
            end
            for dof = 1:15
                % use gpr to set theta_estimates(dof)
                dof
                trial
                activity
                
                meanfunc = @meanZero;
                likfunc = @likGauss;                       
                %meanfunc = @meanConst;
                covfunc = @covSEiso; 
                hyp = struct('mean',[], 'cov', [0,0], 'lik', -1);                            
                % theta_targets(dof)
                hyp2 = minimize(hyp, @gp, -20, @infGaussLik, meanfunc, covfunc, likfunc, v, theta_targets(:, dof));

                hyp2
                [predicted_mean, predicted_var] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, v, theta_targets(:, dof));

                theta_estimates(:, dof) = predicted_mean;
                %calculate add this muscle's mean square error to f                      
                %figure out where f update should be in for loop
                %based on what the mean square loss shoudl be
                %calculatings
                %rn it is the mean square loss over all timesteps
                %for all muscles for 3 trials of each activity
                for i = 1:4000
                    f = f + (theta_estimates(i, dof)-theta_targets(i, dof))^2;
                end
                                   
            end
            theta_estimates
            theta_targets    
        end
    end


    
    f = f/N; 
    f
    A1
    b1
    b2
    alpha1
    din
 
end

function [x,fval,exitflag] = fullexample

%x will be [Ain din y1in y2in]
x0 = [-2,-1,-2,-1,-2,-1,-1,-1, 10, 0.5, 0.5];
%x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
lb = [-3; -3; -3; -3; -3; -3; -3; -3; 0;  -0.99999999999; -0.99999999999];
ub = [0; 0; 0; 0; 0; 0; 0; 0;  20000; 0.99999999999; 0.99999999999];
Aeq = [];
beq = [];
A = [];

b = [];

options = optimoptions('fmincon','Algorithm','sqp');
options.StepTolerance = 1;
options.OptimalityTolerance = 1;
options = optimoptions('fmincon', 'Display', 'iter');

[x,fval,exitflag] = fmincon(@leastSquares,x0,A,b,Aeq,beq,lb,ub,[],options);

end

% function [c, ceq] = confuneq(x)
% c = 0;
% ceq = x(13)*x(14) - x(10);
% 
% end


function expectedJointAngles = calculateExpectedJointAngles(timestep, kinematicsData)
expectedJointAngles = zeros(15, 1);
% kinematicsData is a 4000x69 array (4000 timesteps, 23 markers times 3 xyz
timestepData = kinematicsData(timestep,:); %should be 1x69
[vec1, vec2] = findIndices(23, 21, 17, timestepData);
expectedJointAngles(1) = findAngleBetweenVectors(vec1, vec2); %thumb CMC

[vec1, vec2] = findIndices(21, 17, 18, timestepData);
expectedJointAngles(2) = findAngleBetweenVectors(vec1, vec2); %thumb MCP

[vec1, vec2] = findIndices(17, 18, 19, timestepData);
expectedJointAngles(3) = findAngleBetweenVectors(vec1, vec2); %thumb IP

[vec1, vec2] = findIndices(20, 1, 5, timestepData);
expectedJointAngles(4) = findAngleBetweenVectors(vec1, vec2); %index MCP

[vec1, vec2] = findIndices(1, 5, 6, timestepData);
expectedJointAngles(5) = findAngleBetweenVectors(vec1, vec2); %index PIP

[vec1, vec2] = findIndices(5, 6, 7, timestepData);
expectedJointAngles(6) = findAngleBetweenVectors(vec1, vec2); %index DIP

[vec1, vec2] = findIndices(20, 2, 8, timestepData);
expectedJointAngles(7) = findAngleBetweenVectors(vec1, vec2); %middle MCP

[vec1, vec2] = findIndices(2, 8, 9, timestepData);
expectedJointAngles(8) = findAngleBetweenVectors(vec1, vec2); %middle PIP

[vec1, vec2] = findIndices(8, 9, 10, timestepData);
expectedJointAngles(9) = findAngleBetweenVectors(vec1, vec2); %middle DIP

[vec1, vec2] = findIndices(20, 3, 11, timestepData);
expectedJointAngles(10) = findAngleBetweenVectors(vec1, vec2); %ring MCP

[vec1, vec2] = findIndices(3, 11, 12, timestepData);
expectedJointAngles(11) = findAngleBetweenVectors(vec1, vec2); %ring PIP

[vec1, vec2] = findIndices(11, 12, 13, timestepData);
expectedJointAngles(12) = findAngleBetweenVectors(vec1, vec2); %ring DIP

[vec1, vec2] = findIndices(20, 4, 14, timestepData);
expectedJointAngles(13) = findAngleBetweenVectors(vec1, vec2); %little MCP

[vec1, vec2] = findIndices(4, 14, 15, timestepData);
expectedJointAngles(14) = findAngleBetweenVectors(vec1, vec2); %little PIP

[vec1, vec2] = findIndices(14, 15, 16, timestepData);
expectedJointAngles(15) = findAngleBetweenVectors(vec1, vec2); %little DIP

end

function [v1, v2] = findIndices(marker1, marker2, marker3, timestepData)
% this function finds marker1-marker2 and marker2-marker3, and returns the
% 6 numbers returned in a vector
m1 = marker1*3;
m2 = marker2*3;
m3 = marker3*3;

x1 = timestepData(m1-2) - timestepData(m2-2);
y1 = timestepData(m1-1) - timestepData(m2-1);
z1 = timestepData(m1) - timestepData(m2);
x2 = timestepData(m2-2) - timestepData(m3-2);
y2 = timestepData(m2-1) - timestepData(m3-1);
z2 = timestepData(m2) - timestepData(m3);

v1 = [x1 y1 z1];
v2 = [x2 y2 z2];


end


function angle = findAngleBetweenVectors(v1, v2)

%v1 = [vectors(1) vectors(2) vectors(3)];
%v2 = [vectors(4) vectors(5) vectors(6)];

dotVal = dot(v1, v2);
%dotVal
normVal = norm(v1)*norm(v2);
%normVal

minVal = min(dotVal/normVal, 1);
theta = max(minVal, -1);
angle = real(acosd(theta));

end



