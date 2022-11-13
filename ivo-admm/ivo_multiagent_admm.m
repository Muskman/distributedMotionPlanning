%% Multiagent IVO ADMM: position exchange simulation

% clear the workspace
clear
clc

% environment parameters
n_agents = 4;               % number of agents
n_obstacles = 1;            % number of agents
radius_env = 5;             % size of circular environment 
agents_radius = 0.5;        % radious of agents
obstacles_pose = [0 0];
obstacels_vel = [0 0];
obstacles_radius = 1;

% create an instance of environment
env = env(n_agents,n_obstacles,radius_env,agents_radius,obstacles_pose, obstacels_vel, obstacles_radius);
% env.agents_pose(1,2) = -1; % env.agents_pose(2,2) = -1;
% create instance of agents in the given environment
for index=1:n_agents
    agents(index) = agent(env,index);
end

agents_vel = ones(1,n_agents);
while (norm(agents_vel) > 1e-18)
    while true
        for index = 1:env.n_agents
            [agents(index), env] = agents(index).z_update(env);
            [agents(index), env] = agents(index).x_update(env);
            [agents(index), env] = agents(index).y_update(env);
        end
        [converged, dist] = env.isConverged;
        if converged
            fprintf('ADMM Converged! Convergence Distance %f\n',dist)
            for index = 1:env.n_agents
                agents(index).reinit;
            end
            env.reinit;
            break
        end
    end

    if norm(agents(1).pose - agents(2).pose) < 1e-1
        % keyboard
    end

    for index = 1:env.n_agents
        [env, agents(index)] = agents(index).move(env);
        agents_vel(index) = norm(agents(index).vel);
    end
    env.plot_env();
end