source "./help/geometry_helpers_2d.m"
source "./help/utility.m"
source "./graph_slam.m"

file = "./slam2D_bearing_only_initial_guess.g2o";
file_ground_truth = "./slam2D_bearing_only_ground_truth.g2o";

[num_poses,num_landmarks,id_pose,associations_p_l,associations_p_p,Z_bearing,XR_guess,XL_guess,Zij,id_pose_toindex,rs,ps,land_id,land_pos,landmark_list] = create_initial_guesses(file);

[XR_ground_truth,XL_ground_truth] = load_ground_truth(file_ground_truth);


%Change or remove
num_iterations = 40;
damping = 0.0001;
kernel_threshold = 100;


[XR, XL, chi_stats,chi_stats_p_l,chi_stats_p_p, num_inliers]=LeastSquare(XR_guess, XL_guess, Z_bearing, Zij,
							associations_p_l,
							associations_p_p, landmark_list,
							num_poses, 
							num_landmarks, id_pose,id_pose_toindex,
							num_iterations, 
							damping, 
							kernel_threshold);
              
 %i= plotState(XL, XL_guess, XL_ground_truth)
  
 i= plotStateRobot(XR, XR_guess, XR_ground_truth)


							

