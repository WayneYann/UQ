function  [x_list, fval_list] = surrogate_opt( case_origin, X0 )

%% Optimization
    disp( 'Optimizing the surrogate subspace' );
    rng default % for reproducibility
    %r = generate_sample(1, 1000, 2);
    
    ConstraintFunction = [];
    nvars = case_origin.para.nrxn;
    
    UB = zeros(1,nvars);
    UB(case_origin.para.ActiveId) = 1;
    LB =  -UB;
    options = gaoptimset('Display','iter' ,'UseParallel',true);
    ObjectiveFunction = @(x)multiobjective(x, case_origin);
    options.TolFun = 1.e-4; %default 1.e-4
    options.PopulationSize = 600; %default 200
    options.ParetoFraction = 0.35;
    x_list = {};
    fval_list = {};
    if case_origin.bool.multiobject
        disp( 'GA step #1, limited to active para' );
        options.Generations = 500;
        [x,fval,exitFlag,Output,Population,Score] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
        optPrint( x, Output);
        x_list{end+1} = x;
        fval_list{end+1} = fval;
        
        disp( 'GA step 2, relax the active para' );
        options.Generations = 500;
        UB(:) = 1;
        LB(:) = -1;
        % options.HybridFcn = @fgoalattain;
        % Provide initial population and scores
        options.InitialPopulation = Population;
        options.InitialScores = Score;
        [x,fval,exitFlag,Output,Population,Score] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
        optPrint( x, Output);
        x_list{end+1} = x;
        fval_list{end+1} = fval;
        
        disp( 'GA step 3, decrease TolFun' );
        options.Generations = 5000;
        options.TolFun = 1.e-5;
        % Provide initial population and scores
        options.InitialPopulation = Population;
        options.InitialScores = Score;
        [x,fval,exitFlag,Output,Population,Score] = gamultiobj( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
        optPrint( x, Output);
        x_list{end+1} = x;
        fval_list{end+1} = fval;
    else
        options.Generations = 500;
        options.TolFun = 1.e-5;
        UB(:) = 1;
        LB(:) = -1;
        [x,fval,exitFlag,Output,Population,Score] = ga( ObjectiveFunction, nvars, [], [], [], [], LB, UB, ConstraintFunction, options );
        x_list{end+1} = x;
        fval_list{end+1} = fval;
    end
end

function optPrint( x, Output)
        fprintf('The number of points on the Pareto front was: %d\n', size(x,1));
        fprintf('The average distance measure of the solutions on the Pareto front was: %g\n', Output.averagedistance);
        fprintf('The spread measure of the Pareto front was: %g\n', Output.spread);
end