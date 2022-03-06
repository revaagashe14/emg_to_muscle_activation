

tic
fullExample
toc
t=toc;
    
function vArray = findMuscleActivation(input, startA, endA, startT, endT)
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
    
% 
    b1 = y1in + y2in;
    b2 = y1in * y2in;
    alpha1 = 1 + b1 + b2;
    
    
    vArray = [];
    %din = din/5;
    dataMatrix = load('Data/s1.mat', 'dsfilt_emg', 'finger_kinematics'); 
    %dataMatrix
    for activity = startA:endA 
     
        for trial = startT:endT
                      
            emg_signal = dataMatrix.('dsfilt_emg');
            trial_activity_data = emg_signal{trial, activity};

            u = zeros(4000, 8); % neural activation (used to calc v)
            v = zeros(4000, 8); % muscle activation
                     
            for t = 1:4000  %t is the timestep (out of 4000)              
 
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
    
                        % fval = funcEx(t-din);

                        if (t*5-din >=3) && (t*5-din <= 20000)
                            u(t,x) = alpha1*e(x, round((t*5-din)/5));
                            
                        end
                        if (t-1 >=1)
                            u(t, x) = u(t,x) - b1*u(t-1, x);
                        end

                        if (t-2 >=1)
                            u(t, x) = u(t,x) - b2 * u(t-2, x);
                        end  

                        % calculate v(t,x)
                        v(t,x) = (exp(Ain(x)* u(t,x)) - 1) / (exp(Ain(x))-1);    
                    end % ending muscle for loop
                    
           
             end % ending timestep for loop
             vArray = vertcat(vArray, v);
             
             fprintf("ended timestep loop \n");
        end  % ending activity for loop
    end % ending trial for loop
 
end

function thetaTargetsArray = getThetaTargets(startA, endA, startT, endT)
    thetaTargetsArray = [];
    dataMatrix = load('Data/s1.mat', 'dsfilt_emg', 'finger_kinematics'); 

    for activity = startA:endA 
        for trial = startT:endT
            kinematics = dataMatrix.('finger_kinematics');
            kinematics_data = kinematics{trial, activity};

            theta_targets = zeros(4000, 69);
  
            for t = 1:4000  %t is the timestep (out of 4000)              
                    %calculate theta_targets = correct joint angles for a
                    %given joint positions for a specific trial and activity and timestep
    %                theta_targets(t,:) = calculateExpectedJointAngles(t, kinematics_data);
                     theta_targets(t,:) = kinematics_data(t,:);
            end
            thetaTargetsArray = vertcat(thetaTargetsArray, theta_targets);
        end
    end
end

function GPRmodels = findGPRModels(vArray, thetaTargetsArray)
    GPRmodels = containers.Map('KeyType', 'double', 'ValueType', 'any');
    for dof = 1:69
        dof
        model = fitrgp(vArray, thetaTargetsArray(:, dof), 'Standardize', 1, 'SigmaLowerBound', 5, 'Sigma', 5);
        GPRmodels(dof)=model;
    end
    
end

function errorsEstimatesDictionary = testingFunction(inputs)
   
    trainingStartTrial = 1; % should be 1 --> 5 (or any stop in between)

    trainingEndTrial = 3;

    trainingStartActivity = 1; % should be 1 --> 7

    trainingEndActivity = 7;

    testingStartTrial = 4; % should be 1 --> 5 (or any stop in between)
    testingEndTrial = 5;
    testingStartActivity = 1; % should be 1 --> 7
    testingEndActivity = 7;

    errors = [];
    theta_estimates_total = [];
    
    
    vArrayTraining = findMuscleActivation(inputs,trainingStartActivity, trainingEndActivity, trainingStartTrial, trainingEndTrial);
        
    vArrayTesting = findMuscleActivation(inputs, testingStartActivity, testingEndActivity, testingStartTrial, testingEndTrial);

    size(vArrayTraining)
    size(vArrayTesting)
    
    thetaTargetsTraining = getThetaTargets(trainingStartActivity, trainingEndActivity, trainingStartTrial, trainingEndTrial);
    thetaTargetsTesting = getThetaTargets(testingStartActivity, testingEndActivity, testingStartTrial, testingEndTrial);

    GPRmodels = findGPRModels(vArrayTraining, thetaTargetsTraining);
    
    fprintf("got all gpr models\n");
    
    for dof = 1:69
        modelDof = GPRmodels(dof);
        theta_estimates = predict(modelDof, vArrayTesting); 
        theta_targets_dof = thetaTargetsTesting(:, dof);        
        nrmse = goodnessOfFit(theta_estimates, theta_targets_dof, 'NRMSE');
        errors = [errors, nrmse];
        theta_estimates_total = [theta_estimates_total, theta_estimates];
        
        fprintf("dof is");
        dof
    end
    
    errorsEstimatesDictionary = containers.Map;
    errorsEstimatesDictionary('errors') = errors;
    errorsEstimatesDictionary('estimates') = theta_estimates_total;
    
end

function cost = getCost(inputs)
    global incrementer
    incrementer = incrementer+1;
    incrementer
    
    dict = testingFunction(inputs);
    errorArray = dict('errors');
    cost = mean(errorArray);
    
end

function thetaEstimates = fullExample()
    
    x0 = [-2,-1,-2,-1,-2,-1,-1,-1, 10, 0.5, 0.5];
    lb = [-3; -3; -3; -3; -3; -3; -3; -3; 0;  -0.99999999999; -0.99999999999];
    ub = [0; 0; 0; 0; 0; 0; 0; 0;  20000; 0.99999999999; 0.99999999999];
    Aeq = [];
    beq = [];
    A = [];
    b = [];
    
    global incrementer;
    incrementer = 0;
    
    options = optimoptions('fmincon','Algorithm','sqp');
    options.StepTolerance = 1;
    options.OptimalityTolerance = 1;
    options = optimoptions('fmincon', 'Display', 'iter');

    [x,fval,exitflag] = fmincon(@getCost,x0,A,b,Aeq,beq,lb,ub,[],options);
    fprintf("x is: \n")
    x
    
    dict = testingFunction(x);
    theta_estimates = dict('estimates');
    save('all_estimates.mat', 'theta_estimates');

end