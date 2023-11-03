
clear all;
close all;

theta = nan;
noiseTCorr = nan;

c = questdlg('Run simulated trajectories or open illustrative interface', 'What program should I run?', 'Interface', 'Run trajectories', 'Cancel', 'Interface');



   switch c,
     case 'Interface',
         matlab_visualization
     case 'Run trajectories',
         
         chooseCondition = menu('Choose condition',...
             'Condition 1: theta = 0; noiseTCorr = 20', ...
             'Condition 2: theta = pi/2; noiseTCorr = 20',... 
             'Condition 3: theta = 0; noiseTCorr = 100',...
             'Condition 4: theta = pi/2; noiseTCorr = 100', ...
             'Same as condition 1, but three individuals', ...
             'Same as condition 2, but three individuals');
         
        switch chooseCondition
            case 1
                theta = 0;
                noiseTCorr = 20;
                nIndividuals = 2;
            case 2
                theta = pi/2;
                noiseTCorr = 20;
                nIndividuals = 2;
            case 3
                theta = 0;
                noiseTCorr = 100;
                nIndividuals = 2;
            case 4
                theta = pi/2;
                noiseTCorr = 100;
                nIndividuals = 2;
            case 5
                theta = 0;
                noiseTCorr = 20;
                nIndividuals = 3;
            case 6
                theta = pi/2;
                noiseTCorr = 20;
                nIndividuals = 3;
        end
      apparent_rules_of_a_particle_at_r_theta_from_neighbour(theta, noiseTCorr, nIndividuals);
       case 'Cancel'
   end % switch
   
   
   