source "./help/geometry_helpers_2d.m"

%(minimal) size of pose and landmarks
global pose_dim=3;
global landmark_dim=2;



function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*pose_dim;
endfunction;


function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>199) %num_landmarks
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;



function [e,Jr,Jl]=errorAndJacobian(Xr,Xl,z)
   R=Xr(1:2,1:2);
   t=-R*Xr(1:2,3);
   p_hat = R*Xl +t; 
   z_hat= atan2(p_hat(2), p_hat(1));
   e=z_hat-z;
   

   Jr = zeros(1,3);
  Jl = zeros(1,2);
  
  x_hat = p_hat(1); y_hat = p_hat(2);
  J_atan = 1/(x_hat^2 + y_hat^2)*[-y_hat x_hat];
  Jr = J_atan*[-R R*[Xl(2) -Xl(1)]'];
  Jl = J_atan*R;
endfunction;



function [eij,Ji,Jj]=errorAndJacobianPosePose(Xi,Xj,Z)
   Ri=Xi(1:2,1:2);
   ti=Xi(1:2,3);
   Rj=Xj(1:2,1:2);
   tj=Xj(1:2,3);

  tij=tj-ti;
  Ri_transpose=Ri';
  Ji=zeros(6,3);
  Jj=zeros(6,3);
  
  dR_0 = [0 -1;1 0];
  dR = Ri_transpose*dR_0*Rj;
  
  Jj(1:4,3)=reshape(dR, 4, 1);
  Jj(5:6,1:2)=Ri_transpose;
  
  Jj(5:6,3)=Ri_transpose*dR_0*tj;
  Ji=-Jj;

  Z_hat=eye(3);
  Z_hat(1:2,1:2)=Ri_transpose*Rj; 
  Z_hat(1:2,3)=Ri_transpose*tij;
  eij=flattenIsometryByColumns(Z_hat-Z);
   
endfunction;



function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl;
  endfor;
endfunction;



function [XR, XL, chi_stats,chi_stats_p_l,chi_stats_p_p, num_inliers]=LeastSquare(XR, XL, Z, Zij, 
							associations_p_l,
              associations_p_p, landmark_list,
							num_poses, 
							num_landmarks,id_pose, id_pose_toindex,
							num_iterations, 
							damping, 
							kernel_threshold
              )
  global pose_dim;
  global landmark_dim;
  % Information Matrices for the Z and Zij measurement respectively
  sigma_bearing = 0.0002;
  sigma_p_p = eye(6)*6;
  poses=[];
  chi_stats=zeros(1,num_iterations);
  chi_stats_p_p = zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  % size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*199; 
  for (iteration=1:num_iterations)
  
      printf("iteration %d\n",iteration)  
      
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
    H_p=zeros(system_size, system_size);
    b_p=zeros(system_size,1);
        H_l=zeros(system_size, system_size);
    b_l=zeros(system_size,1);
    chi_stats(iteration)=0;
    % First compute the error and accumulate in H and b all the terms due the bearing measurement error 
    for (measurement_num=1:size(Z,2))
      pose_id=associations_p_l(1,measurement_num);
      landmark_index=associations_p_l(2,measurement_num);
  if length(find(landmark_list==landmark_index))==0 || XL(:,landmark_index)==[0,0]
    continue;
  endif
      poses(end+1)=pose_id;
        z=Z(:,measurement_num);
      Xr=XR(:,:,id_pose_toindex(pose_id));
      Xl=XL(:,landmark_index);
      [e,Jr,Jl] = errorAndJacobian(Xr, Xl, z);
      
      
   
      gamma_ = 1;


        chi=e'*e;
        if (chi>kernel_threshold)
          e*=sqrt(kernel_threshold/chi);
          chi=e'*e;
          %chi=kernel_threshold;
        else
          num_inliers(iteration)++;
        end;
         
      chi_stats(iteration)+=chi;

      Hrr=Jr'*gamma_*sigma_bearing*Jr;
      Hrl=Jr'*gamma_*sigma_bearing*Jl;
      Hll=Jl'*gamma_*sigma_bearing*Jl;
      br=Jr'*gamma_*sigma_bearing*e;
      bl=Jl'*gamma_*sigma_bearing*e;



      pose_matrix_index=poseMatrixIndex(id_pose_toindex(pose_id), num_poses, num_landmarks);
      landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

      H_l(pose_matrix_index:pose_matrix_index+pose_dim-1,
	pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;

      H_l(pose_matrix_index:pose_matrix_index+pose_dim-1,
	landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;

      H_l(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
	landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;

      H_l(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
	pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrl';

      b_l(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
      b_l(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;

    endfor

    chi_stats_p_l(iteration) = chi_stats(iteration);

  %POSES
    % Then compute the error and keep accumulating in H and b all the terms due the relative poses measurement error
    for (measurement_num=1:size(Zij,3))
 
      pose_index_i=associations_p_p(1,measurement_num);
      pose_index_j=associations_p_p(2,measurement_num);
       % if length(find(poses==pose_index_i))==0
    %continue;
  %endif
      z=Zij(:,:,measurement_num);
      Xri=XR(:,:,id_pose_toindex(pose_index_i));
      Xrj=XR(:,:,id_pose_toindex(pose_index_j));
      [e,Jri,Jrj] = errorAndJacobianPosePose(Xri, Xrj, z);
      gamma_ = 1;
  
          chi=e'*e;
        if (chi>kernel_threshold)
          e*=sqrt(kernel_threshold/chi);
             chi=e'*e;
          %chi=kernel_threshold;
        else
          num_inliers(iteration)++;
        end;
 
      chi_stats(iteration)+=chi;
      chi_stats_p_p(iteration)+=chi;

      Hii=Jri'*gamma_*sigma_p_p*Jri;
      Hij=Jri'*gamma_*sigma_p_p*Jrj;
      Hjj=Jrj'*gamma_*sigma_p_p*Jrj;
      bri=Jri'*gamma_*sigma_p_p*e;
      brj=Jrj'*gamma_*sigma_p_p*e;



      pose_matrix_index_i=poseMatrixIndex(id_pose_toindex(pose_index_i), num_poses, num_landmarks);
      pose_matrix_index_j=poseMatrixIndex(id_pose_toindex(pose_index_j), num_poses, num_landmarks);

      H_p(pose_matrix_index_i:pose_matrix_index_i+pose_dim-1,
  pose_matrix_index_i:pose_matrix_index_i+pose_dim-1)+=Hii;

      H_p(pose_matrix_index_i:pose_matrix_index_i+pose_dim-1,
  pose_matrix_index_j:pose_matrix_index_j+pose_dim-1)+=Hij;

      H_p(pose_matrix_index_j:pose_matrix_index_j+pose_dim-1,
  pose_matrix_index_j:pose_matrix_index_j+pose_dim-1)+=Hjj;

      H_p(pose_matrix_index_j:pose_matrix_index_j+pose_dim-1,
  pose_matrix_index_i:pose_matrix_index_i+pose_dim-1)+=Hij';

      b_p(pose_matrix_index_i:pose_matrix_index_i+pose_dim-1)+=bri;
      b_p(pose_matrix_index_j:pose_matrix_index_j+pose_dim-1)+=brj;

    endfor
    H=H_p;
    b=b_p;
      if (num_landmarks) 
       H+=H_l;
       b+=b_l;
    endif;
  

    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);


    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  endfor
endfunction


function v=flattenIsometryByColumns(T)
v=zeros(6,1);
v(1:4)=reshape(T(1:2,1:2),4,1);
v(5:6)=T(1:2,3);
endfunction

% plot landmarks and poses
%
%
%
function i = plotState(XL, XL_guess, XL_gt)
%plot landmarks
hold on;
plot(XL(1,:),XL(2,:),'b*',"linewidth",2);
hold on;
plot(XL_guess(1,:),XL_guess(2,:),'ro',"linewidth",2);
hold on;
plot(XL_gt(1,:),XL_gt(2,:),'g*',"linewidth",2);
hold on;
legend("estimate","initial guess","ground truth")
i = 1;
endfunction

function i = plotStateRobot(XR, XR_guess, XR_gt)
%plot landmarks
hold on;
plot(XR(1,3,:),XR(2,3,:),'r+',"linewidth",2);
hold on;
plot(XR_guess(1,3,:),XR_guess(2,3,:),'bo',"linewidth",2);
hold on;
plot(XR_gt(1,3,:),XR_gt(2,3,:),'g',"linewidth",2);
hold on;
legend("estimate","initial guess","ground truth")
i = 1;
endfunction
