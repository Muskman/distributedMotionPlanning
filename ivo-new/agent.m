classdef agent 

    properties
        index               % agent index
        radius 
        pose
        vel
        vel_max

        goal
        r_obstacles
        pose_obstacles
        vel_obstacles
        safety_margin       % collision safety margin

        u                   % control

        lambda              % smoothness parameter
        avoidance_range     % obstacle avoidance range
    end

    methods 
        function default_agent = agent(env, index)
            % creats a default agent based on the 
            default_agent.index = index;
            default_agent.radius = env.agents_radius;
            default_agent.pose = env.agents_pose(index,:);
            default_agent.vel = [0 0];
            default_agent.vel_max = 2;

            default_agent.goal = env.agents_goal(index,:);
            default_agent.r_obstacles = [repmat(env.agents_radius, env.n_agents, 1); repmat(env.obstacles_radius, env.n_obstacles, 1)]; % add radius of static obstacles that are in the env
            
            pose_obstacles = env.agents_pose;
            pose_obstacles(index,:) = [];
            default_agent.pose_obstacles = [pose_obstacles; env.obstacles_pose]; % add pose of static obstacles that are in the env

            vel_obstacles = env.agents_vel;
            vel_obstacles(index,:) = [];
            default_agent.vel_obstacles = [vel_obstacles; zeros(size(env.obstacles_pose))]; % add vel of static obstacles that are in the env

            default_agent.safety_margin = 2*default_agent.radius;

            default_agent.u = (default_agent.goal-default_agent.pose)*eps;
            default_agent.lambda = 4;
            default_agent.avoidance_range = 6;
        end

        function [c,ceq] = generate_constraints(agent,u)
            % generate constraints
            n_obs = size(agent.pose_obstacles,1);
            c = zeros(1,n_obs);
            for i = 1:n_obs
                p_relative = (agent.pose_obstacles(i,:)-agent.pose); v_relative = (agent.vel_obstacles(i,:)-agent.vel); 
                if p_relative*v_relative' <= norm(p_relative)*norm(v_relative)*cosd(80) && norm(p_relative) <= agent.avoidance_range
                    critical_dist = norm(p_relative)^2 - (agent.radius + agent.safety_margin + agent.r_obstacles(i))^2;
                    c(i) = ((p_relative*(v_relative-u)')/norm(v_relative-u))^2 - critical_dist;
                else
                    c(i) = 0;
                end
                ceq = [];
            end
        end

        function [new_env, new_agent] = move_and_update(agent,env)
            
            p_obstacles = env.agents_pose;
            p_obstacles(agent.index,:) = [];
            agent.pose_obstacles(1:env.n_agents-1,:) = p_obstacles;

            v_obstacles = env.agents_vel;
            v_obstacles(agent.index,:) = [];
            agent.vel_obstacles(1:env.n_agents-1,:) = v_obstacles;
            
            if norm(agent.pose-agent.goal) > 1e-1
                vel_desired = agent.vel_max*(agent.goal-agent.pose)/norm(agent.goal-agent.pose);
                objective = @(u) (norm(vel_desired-(agent.vel+u)))^2 + agent.lambda * (norm(u)^2);
                constraints = @(u) agent.generate_constraints(u);

                problem.objective = objective;
                problem.nonlcon = constraints;

                options = optimoptions('fmincon','display','off','Algorithm','sqp');
                problem.options = options;
                problem.solver = 'fmincon';
                problem.x0 = agent.u;
                
                agent.u = fmincon(problem);
                
                v = agent.vel + agent.u; 
                if norm(v)>agent.vel_max
                    agent.vel = agent.vel_max*v/norm(v);
                else
                    agent.vel = v;
                end
            else
                agent.vel = 0;
                agent.u = 0;
            end

            agent.pose = agent.pose + agent.vel*env.deltaT;

            new_env = env;
            new_env.agents_pose(agent.index,:) = agent.pose;
            new_env.agents_vel(agent.index,:) = agent.vel;

            new_agent = agent;
        end

    end

end