%leastSquares(0,0,0,0,0,0,0,0,0,0,0)
fullexample
function f = leastSquares(A1)
    
    
    A2 = A1(2);
    A3 = A1(3);
    A4 = A1(4);
    A5 = A1(5);
    A6 = A1(6);
    A7 = A1(7);
    A8 = A1(8);
    din = A1(9);
    y1in = A1(10);
    y2in = A1(11);
    A1 = A1(1);
    Ain = [A1; A2; A3; A4; A5; A6; A7; A8];
  
    
    
    t_end = 4000;

    b1 = y1in + y2in;
    b2 = y1in * y2in;
    alpha1 = 1 + b1 + b2;
    
    u = zeros(4000, 8);
    v = zeros(4000, 8);
   
    f = 0;
    
    dataMatrix = load('Data/s1.mat', 'dsfilt_emg', 'finger_kinematics'); 
    %dataMatrix
    for trial = 1:3 %should be 3
       % trial
        for activity = 1:7 %should be 7
            %activity
                for t = 1:t_end  %t is the timestep (out of 4000)
                    

                    emg_signal = dataMatrix.('dsfilt_emg');
                    trial_activity_data = emg_signal{trial, activity};
                    
                    kinematics = dataMatrix.('finger_kinematics');
                    kinematics_data = kinematics{trial, activity};
                    
                    %calculate theta_targets
                    theta_targets = calculateExpectedJointAngles(t, kinematics_data);
                    theta_targets
                    
                    N = length(theta_targets); %should be 8 for 8 muscles for the orig dataset
                    theta_estimates = zeros(N);  %use gpr to find this
                   % trial_activity_data

                % calculating e below
                    for x = 1:8  %x is the muscle (out of 8)
                        %x
                        e(t, x) = trial_activity_data(t, x);  % change this to get it from data
                        
                        
                        
                        %evaluye = e(t,x)
                        
                        
                        % find u(t, x): (later account for base cases: what if din>t? are the if cases right)
                        if (t-din >=1) && (t-din <= t)
                            u(t, x) = u(t,x) + alpha1*e(round(t-din), x);
                        end
                        if (t-1 >=1)
                            u(t, x) = u(t,x) - b1*u(t-1, x);
                        end

                        if (t-2 >=1)
                            u(t, x) = u(t,x) - b2 * u(t-2, x);
                        end             


                        % calculate v(t,x)
                        v(t,x) = (exp(Ain(x)* u(t,x)) - 1) / (exp(Ain(x))-1);  %why is htis NaN
                        a = v(t,x);
                        
                        b = Ain(x);
                        
                        
                        c = u(t,x);
                        
                    end
                    
                    for dof = 1:15
                        % use gpr to set theta_estimates(dof)
                        meanfunc = @meanZero;
                        likfunc = @likGauss;
                        
                        %meanfunc = @meanConst;
                        covfunc = @covSEiso; 

                        hyp = struct('mean', [], 'cov', [0,0], 'lik', 0);
                            
                        
                        %theta_targets(dof)
                        hyp2 = minimize(hyp, @gp, -20, @infGaussLik, meanfunc, covfunc, likfunc, v(t), theta_targets(dof));

                       
                        [predicted_mean, predicted_var] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, v(t), theta_targets(dof));

                        theta_estimates(dof) = predicted_mean;
                        predicted_mean

                        %calculate add this muscle's mean square error to f
                        
                        
                        %figure out where f update should be in for loop
                        %based on what the mean square loss shoudl be
                        %calculating
                        %rn it is the mean square loss over all timesteps
                        %for all muscles for 3 trials of each activity
                     
                        f = f + (theta_estimates(x)-theta_targets(x))^2;

                    end

                end
        end
    end
    

    f = f/N; 
 
end

function [x,fval,exitflag] = fullexample

%x will be [Ain din y1in y2in]
x0 = [-2,-1,-2,-1,-2,-1,-1.567,-1, 3, -0.123,0.75];
%x0 = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
lb = [-3; -3; -3; -3; -3; -3; -3; -3; 0;  -0.99999999999; -0.99999999999];
ub = [0; 0; 0; 0; 0; 0; 0; 0;  Inf; 0.99999999999; 0.99999999999];
Aeq = [0 0 0 0 0 0 0 0 0 0 0 ];
beq = 0;
A = [0 0 0 0 0 0 0 0 0 0 0  
     0 0 0 0 0 0 0 0 0 0 0 
     0 0 0 0 0 0 0 0 0 0 0 ];
b = [0; 0; 0];

options = optimoptions('fmincon','Algorithm','active-set');

[x,fval,exitflag] = fmincon(@leastSquares,x0,A,b,Aeq,beq,lb,ub, [], options);

end

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



