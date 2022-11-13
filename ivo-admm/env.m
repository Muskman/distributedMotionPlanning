classdef env
    properties
        deltaT          % time scale resolution
        t               % current simulation time counter
        T               % time horizon
        
        radius_env      % radius of environment
        n_agents        % number of agets
        n_obstacles     % number of static obstacles

        agents_radius   % agents radius
        agents_pose     % agents cartesian positions
        agents_vel      % agents velocity
        agents_goal     % agents goal positions
        agents_colors   % agent colors
        agents_trajX    % agent trajectories x coordinate
        agents_trajY    % agent trajectories y coordinate

        agents_xN       % agents estimates of other agents control
        agents_yN       % agents estimates of other agents dual variable
        agents_zN       % agents control iterate available for sharing
        
        agents_u        % agents control to be executed
        
        obstacles_pose  % obstacle pose (currently static)
        obstacels_vel   % obstacle velocity
        obstacles_radius% obstacle radius (same for all obstacles)
    end

    methods
        function default_env = env(n_agents,n_obstacles,radius_env,agents_radius,obstacles_pose, obstacels_vel,obstacles_radius)
            default_env.deltaT = 0.1;
            default_env.t = 1;
            default_env.T = 1000;
            
            default_env.radius_env = radius_env;
            default_env.n_agents = n_agents;
            default_env.n_obstacles = n_obstacles;
            
            default_env.agents_radius = agents_radius;

            agents_angle = 0:2*pi/n_agents:2*pi*(1-1/n_agents); % agent angular positions
            default_env.agents_pose = radius_env*[cos(agents_angle') sin(agents_angle')];

            default_env.agents_vel = zeros(n_agents,2);

            default_env.agents_goal = zeros(size(default_env.agents_pose)); % goal for exchanging positions
            default_env.agents_goal(1:n_agents/2,:) = default_env.agents_pose((1:n_agents/2)+n_agents/2,:);
            default_env.agents_goal((n_agents/2+1):n_agents,:) = default_env.agents_pose(((n_agents/2+1):n_agents)-n_agents/2,:);

            default_env.agents_colors = jet(n_agents);
            default_env.agents_trajX = repmat(default_env.agents_pose(:,1),1,default_env.T);
            default_env.agents_trajY = repmat(default_env.agents_pose(:,2),1,default_env.T);

            default_env.agents_xN = repmat(zeros(n_agents,2),[1 1 n_agents]);
            default_env.agents_yN = repmat(zeros(n_agents,2),[1 1 n_agents]);
            default_env.agents_zN = zeros(n_agents,2);

            default_env.agents_u = zeros(n_agents,2);    

            default_env.obstacles_pose = obstacles_pose;
            default_env.obstacels_vel = obstacels_vel; 
            default_env.obstacles_radius = obstacles_radius;
        end

        function new_env = update_env(env,agents)
            new_env = env;
            for index = 1:env.n_agents
                new_env.agents_pose(index,:) = agents(index).pose;
                new_env.agents_vel(index,:) = agents(index).vel;
            end
            new_env.t = new_env.t + 1;
            new_env.agents_trajX(:,new_env.t) = new_env.agents_pose(:,1);
            new_env.agents_trajY(:,new_env.t) = new_env.agents_pose(:,2);
        end

        function [converged, dist] = isConverged(env)
            dist = (1/env.n_agents)*vecnorm(vecnorm(vecnorm(env.agents_xN-env.agents_zN,2,2)),2,3);
            fprintf('Convergence Distance %f\n',dist)
            converged = dist < 1e-3;
        end

        function new_env = reinit(env)
            n = env.n_agents;
            new_env = env;

            new_env.agents_xN = repmat(zeros(n,2),[1 1 n]);
            new_env.agents_yN = repmat(zeros(n,2),[1 1 n]);
            new_env.agents_zN = zeros(n,2);
            % new_env.agents_u = zeros(n,2);
        end

        function [] = plot_env(env)
            % plot the current state of the environment
            for i = 1:env.n_agents
                filledCircle(env.agents_pose(i,:),env.agents_radius,100,env.agents_colors(i,:));
                %plot(env.agents_pose(i,1), env.agents_pose(i,2),'o','MarkerEdgeColor','k','MarkerFaceColor',env.agents_colors(i,:),'MarkerSize',env.agents_radius/0.035278)
                hold on; axis square;
                plot(env.agents_trajX(i,1:env.t), env.agents_trajY(i,1:env.t),'-','Color',env.agents_colors(i,:),'LineWidth',2)
            end

            n_obs = size(env.obstacles_pose,1);
            for i = 1:n_obs
                filledCircle(env.obstacles_pose(i,:),env.obstacles_radius,100,'k');
                %plot(env.obstacles_pose(i,1),env.obstacles_pose(i,2),'o','MarkerFaceColor','k','MarkerSize',env.obstacles_radius/0.035278)
            end

            xlim([-env.radius_env*2 env.radius_env*2])
            ylim([-env.radius_env*2 env.radius_env*2])    
            drawnow;
            hold off
        end
    end
end


function h = filledCircle(center,r,N,color)
% FILLEDCIRCLE Filled circle drawing
% 
% filledCircle(CENTER,R,N,COLOR) draws a circle filled with COLOR that 
% has CENTER as its center and R as its radius, by using N points on the 
% periphery.
%
% Usage Examples,
%
% filledCircle([1,3],3,1000,'b'); 
% filledCircle([2,4],2,1000,'r');
    THETA=linspace(0,2*pi,N);
    RHO=ones(1,N)*r;
    [X,Y] = pol2cart(THETA,RHO);
    X=X+center(1);
    Y=Y+center(2);
    h=fill(X,Y,color);
    axis square;
end
