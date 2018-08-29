function  out  = desired_states( arg_steer_angle)

long_speed = 80/3.6;
desired_yaw_rate =  long_speed*arg_steer_angle/(2.5 + (long_speed^2)*0.132);


out = [0, 0, desired_yaw_rate, 0, 0 ];

end
