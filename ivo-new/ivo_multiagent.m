%% Multiagent IVO: position exchange simulation

% clear the workspace
clear
clc

% environment parameters
n_agents = 2;               % number of agents
n_obstacles = 1;            % number of agents
radius_env = 5;             % size of circular environment 
agents_radius = 0.5;        % radious of agents
obstacles_pose = [0 0];
obstacles_radius = 1;

% create an instance of environment
env = env(n_agents,n_obstacles,radius_env,agents_radius,obstacles_pose,obstacles_radius);

% create instance of agents in the given environment
for index=1:n_agents
    agents(index) = agent(env,index);
end

agents_vel = ones(1,n_agents);
while (norm(agents_vel) > 1e-8)
    for index = 1:env.n_agents
        [~, agents(index)] = agents(index).move_and_update(env);
        agents_vel(index) = norm(agents(index).vel);
    end
    env = env.update_env(agents);
    env.plot_env();
end