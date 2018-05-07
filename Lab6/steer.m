function q_new_coord = steer(q_rand, q_nearest, mini_distance, eps)
   qnew = [0 0];
   
   % Steer towards q_rand with maximum step size of eps
   if mini_distance >= eps
       qnew(1) = q_nearest(1) + ((q_rand(1)-q_nearest(1))*eps)/norm(q_rand - q_nearest);
       qnew(2) = q_nearest(2) + ((q_rand(2)-q_nearest(2))*eps)/norm(q_rand - q_nearest);
   else
       qnew(1) = q_rand(1);
       qnew(2) = q_rand(2);
   end   
   q_new_coord = [qnew(1), qnew(2)];
end