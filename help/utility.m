source "./help/geometry_helpers_2d.m"

% This is the function that parse the file containing the initial guesses on robot poses and odometry measurement and range measurement and build 
% the correct structures in order to launch the LS solver

function [num_poses,num_landmark,id_pose,associations_p_l,associations_p_p,Z_bearing,XR_guess,XL_guess,Zij,id_pose_toindex,rs,ps,land_id,land_pos,landmark_list] = create_initial_guesses(file) 
	f = fopen(file);
	initial_guess_poses = [];
	% vector of association for the bearing measurement
	associations_p_l = [];
	% vector of association for the odometry measurement
	associations_p_p = [];
	Z_range = [];
  Z_bearing=[];
	%Zij will contain the relative displacement between position i and j
	Zij = zeros(3,3,1);


	%bookkeeping: to and from mapping between robot pose (x,y, theta) and state vector and landmark indices (i) and state vector and viceversa
	flag = 0;
  flag2 = 0;
  id_pose_toindex=[];
  id_pose = [];
  landmark_list=[];
	land_meas_counter = [];
	l = fgetl(f);
	while ischar(l)
		d = strsplit(l);
		header = d{1};
		%get all initial guesses for the robot poses
		if strcmp(header,"VERTEX_SE2") == 1 %POSE
			id_pose(end+1) = str2double(d{2});
      id_pose_toindex(str2double(d{2})) = length(id_pose);
			x_r = str2double(d{3});
			y_r = str2double(d{4});
			theta_r = str2double(d{5});
     
      if flag2 == 0
				XR_guess(:,:,1) = v2t([x_r;y_r;theta_r]);
				flag2 = 1;
			else	
				XR_guess(:,:,end+1) = v2t([x_r;y_r;theta_r]);
			end
		end
		if strcmp(header,"EDGE_SE2") == 1 %TRANSITION
			% when we are parsing the odometry measurements then we need to build the virtual measurement between pose j and pose i

			id_pose_i = str2double(d{2});
			id_pose_j = str2double(d{3});
      
			dx_r = str2double(d{4});
			dy_r = str2double(d{5});
			dtheta_r = str2double(d{6});
			dXij = v2t([dx_r;dy_r;dtheta_r]); %TODO: replace with bearing
			associations_p_p(:,end+1) = [id_pose_i;id_pose_j];
			if flag == 0
				Zij(:,:,1) = dXij;
				flag = 1;
			else
				Zij(:,:,end+1) = dXij;
			end
		end
    
    
    if strcmp(header,"EDGE_BEARING_SE2_XY") == 1 %OBSERVATION
			% When I am parsin the measurements I only retain the landmarks for which I have at least 2 range measurements, otherwise I can't compute
			% a correct initialization of the landamrks.
      %i can still switch to extracting all landmarks by just uncommenting the below script
			state_pose = str2double(d{2});
			id_landmark = str2double(d{3});
     if id_landmark >0
      if length(land_meas_counter) < id_landmark 
				% In this case I have never encountered the landmark so far so I initialize the counter to 1
				land_meas_counter(id_landmark) = 1;
        %uncomment to keep track of all seen landmarks (despite being seen only once)
                %landmark_list(end+1) = id_landmark;

			elseif land_meas_counter(id_landmark) == 0
				% In this case i have already encountered a landmark whose id was bigger than current landmark_id and so land_meas_counter(id_landmark) was
				% set to 0 and I need to initialize it to 1
				land_meas_counter(id_landmark) = 1;
                %uncomment to keep track of all seen landmarks (despite being seen only once)

                %landmark_list(end+1) = id_landmark;


			else
				% I have already encountered this landmark and I simply increase the counter
				land_meas_counter(id_landmark) = land_meas_counter(id_landmark) + 1;

       
    end
    
    	if land_meas_counter(id_landmark) == 2
          landmark_list(end+1) = id_landmark;				
 
			
			end
          Z_bearing(end+1) = str2double(d{4});
			    associations_p_l(:,end+1) = [state_pose;id_landmark];

		
      endif
		end
		l = fgetl(f);
	end


	num_poses = length(id_pose);
	num_landmark = length(landmark_list);


	XL_guess = zeros(2,num_landmark);

  
	for (lan_num=1:num_landmark)
		
    [rs,ps,zs, land_id] = find_bearings(landmark_list(lan_num),associations_p_l,Z_bearing,XR_guess, id_pose_toindex,Zij);
		land_pos = triangulate(ps,rs,zs,Zij);
		XL_guess(:,landmark_list(lan_num)) = land_pos;
	end

end

%Extracting the ground truth
%this is a more straight forward approach of the Load2Go library
function [XR_ground_truth,XL_ground_truth] = load_ground_truth(file)
	XL_ground_truth = [];
	XR_ground_truth = [];


	flag = 0;
 
	f = fopen(file);
	l = fgetl(file);
	while ischar(l)
		d = strsplit(l);
		header = d{1};
		%get all initial guesses for the robot poses
		if strcmp(header,"VERTEX_XY") == 1
			id_landmark = str2double(d{2});
			x_l = str2double(d{3});
      
			y_l = str2double(d{4});
      if id_landmark>0
			XL_ground_truth(:,id_landmark) = [x_l;y_l];
      end;
    
		end
		if strcmp(header,"VERTEX_SE2") == 1
			x_r = str2double(d{3});
			y_r = str2double(d{4});
			theta_r = str2double(d{5});
			if flag == 0
				XR_ground_truth(:,:,1) = v2t([x_r;y_r;theta_r]);
				flag = 1;
			else	
				XR_ground_truth(:,:,end+1) = v2t([x_r;y_r;theta_r]);
			end
		end
    


		l = fgetl(f);
	end

	end




%	ps: (2xnum_poses) vector of poses from which the landmark has been observed
%	rs: (1xnum_poses) vector of ranges measurement made by poses i in ps with the correct order
% output
%	l: a guess for the landmark position

%This might require a modification as might not be correct
function l = triangulate(ps,rs,zs,Zij)

  
  
  

  robot_pose = v2t([ps(1,1) ps(2,1) ps(3,1)]); %robot_pose ==ps ..i have this
  
 
    
     u_x = zs(1,1);
     u_theta = zs(3,1);
     
     
 
     

     delta_pose = v2t([u_x 0 u_theta]);  %get Zij of that robot... and calculate position with respect to world ...%do it in find landmark
     robot_pose = robot_pose * delta_pose; %position of the robot w.r.t the word in the second observation (l_2_2)
     z_ob1= rs(1)
     z_ob2=rs(2)


        if abs(z_ob1)<0.0001 || abs(u_x)<0.2 %if the first orientation or u_x is 0 discard the candidate landmark since a triangulation is not possible due to
                                             %the chosen formula
         l=[0,0]
       
        elseif abs(z_ob2-z_ob1)<0.03 % the observation of the landmark is almost on the same line: the triangulation is not possible %just make sure that the diff of orentation is not zero
        l=[0,0]
       else
       a11 = delta_pose(1,1);
       a12 = delta_pose(1,2);
       a21 = delta_pose(2,1);
       a22 = delta_pose(2,2);
       u_x = delta_pose(1,3);

       A = a11*cos(z_ob2) + a12*sin(z_ob2);
       B = a21*cos(z_ob2) + a22*sin(z_ob2);

       k = -u_x*sin(z_ob1) / (sin(z_ob1)*A - cos(z_ob1)*B);
       l_2_2 = [k*cos(z_ob2);k*sin(z_ob2)];

        
        %update triangulated_lm (shift the lm in the word frame)
        lm_world = robot_pose * [l_2_2;1];
        l = lm_world(1:2,1);
        endif
end

% function that finds all the ranges measurement associated to a landmark
% input:
%	land_id: id of the target landmark
% 	ass: 2xnum_measurements. 
%   ass(:,k)=[p_id,l_id]' means the kth measurement
%   refers to an observation made from pose p_id, that
%   observed landmark l_id
% 	b_meas: vector of all the range measurements
% 	XR_guess: (3,3,num_poses) initial guesses of the robot poses
% output:
%	rs: vetcor of range measurements associated to land_id
%	px: vector of the poses from which the measurements in rs were taken


function [triangulated_lm_id, triangulated_lm] = triangulation(poses, transitions, observations)  #consecutive observations needs to have some lm seen
  triangulated_lm_id = []; % the lm id of the triangulated landmark
  triangulated_lm = [];
  robot_pose = v2t([poses(2).x poses(2).y poses(2).theta]); %robot_pose ==ps ..i have this
  
  for i=1:length(observations)-1
     ob1 = observations(i).observation;
     ob2 = observations(i+1).observation; % obesrvations Z_bearings ..in this case rs ..i have this
     trans = transitions(i+1).v;  %transitionts Zij ..i have this
     u_x = trans(1);
     u_theta = trans(3);
     
     
     %find the landmark to triangulate ..i have this
     ids_ob1 = zeros(1,length(ob1));
     ids_ob2 = zeros(1,length(ob2));
     for j=1:length(ob1)
      ids_ob1(j) = ob1(j).id;
     endfor
     for j=1:length(ob2)
      ids_ob2(j) = ob2(j).id;
     endfor
     
     adiacent_ids=intersect(ids_ob1,ids_ob2); %find common landmark ..i have this
     lm_to_triangulate = setdiff(adiacent_ids,triangulated_lm_id);
     
     delta_pose = v2t([u_x 0 u_theta]);  %get Zij of that robot... and calculate position with respect to world ...%do it in find landmark
     robot_pose = robot_pose * delta_pose; %position of the robot w.r.t the word in the second observation (l_2_2)
     
     %compute the landmarks to triangulate
     for j=1:length(lm_to_triangulate)
        lm_id = lm_to_triangulate(j);
        lm = ones(3,1); %HOMOGENEOUS COORDINATES
        %find measurements
        z_ob1_index = find(ids_ob1==lm_id);
        z_ob1 = ob1(z_ob1_index).bearing;
        if abs(z_ob1)<0.0001 || abs(u_x)<0.2 %if the first orientation or u_x is 0 discard the candidate landmark since a triangulation is not possible due to
                                             %the chosen formula
          %printf("discarded landmark %d for these coupled poses\n",lm_id) %just make sure that the orentation is not zero
          continue
        endif
        z_ob2_index = find(ids_ob2==lm_id);
        z_ob2 = ob2(z_ob2_index).bearing;
        if abs(z_ob2-z_ob1)<0.03 % the observation of the landmark is almost on the same line: the triangulation is not possible %just make sure that the diff of orentation is not zero
          continue  
        endif
        a11 = delta_pose(1,1);
        a12 = delta_pose(1,2);
        a21 = delta_pose(2,1);
        a22 = delta_pose(2,2);
        u_x = delta_pose(1,3);

        A = a11*cos(z_ob2) + a12*sin(z_ob2);
        B = a21*cos(z_ob2) + a22*sin(z_ob2);

        k = -u_x*sin(z_ob1) / (sin(z_ob1)*A - cos(z_ob1)*B);
        l_2_2 = [k*cos(z_ob2);k*sin(z_ob2)];

        
        %update triangulated_lm (shift the lm in the word frame)
        lm_world = robot_pose * [l_2_2;1];
        triangulated_lm(:,end+1) = lm_world(1:2,1);
        triangulated_lm_id(end+1) = lm_id;
     endfor
  endfor
endfunction






function [rs,px,zx, land_id] = find_bearings(land_id,ass,b_meas,XR_guess, id_pose_toindex, Zij)
  	rs = [];
	px = [];
	for(i=1:length(b_meas))
		l_id = ass(2,i);

		if l_id == land_id
      disp("land id ");
                disp(land_id);

			z = b_meas(i);
			rs(end+1) = z;
			pose_id = ass(1,i);
			Xr = XR_guess(:,:,id_pose_toindex(pose_id));
      if i>1
      Zr= Zij(:,:,id_pose_toindex(pose_id)-1);
        euclidian_state2 = t2v(Zr); 
			zx(:,end+1) = euclidian_state2(1:3);
      end
			euclidian_state = t2v(Xr); 
			px(:,end+1) = euclidian_state(1:3);
    
		end
	end
end