subplot(2,2,1)
tempo = steer(:,1);
plot( tempo, x(1,:)');
hold on
plot( tempo, ones(size(tempo))*pi/18)
plot( tempo, -ones(size(tempo))*pi/18)
title('Deslizamento lateral')
xlabel('tempo (s)')
ylabel('\beta (rad)')

subplot(2,2,2)
plot( steer(:,1), x(2,:)');
hold on
plot( steer(:,1), desired_yaw_rate);
plot( tempo, ones(size(tempo))*pi/6)
plot( tempo, -ones(size(tempo))*pi/6)
title('Taxa de guinada')
xlabel('tempo (s)')
ylabel('d\psi/dt (rad/s)')

subplot(2,2,3)
plot( steer(:,1), x(4,:)');
hold on
plot( tempo, ones(size(tempo))*pi/18)
plot( tempo, -ones(size(tempo))*pi/18)
title('Rolagem')
xlabel('tempo (s)')
ylabel('\phi (rad)')

subplot(2,2,4)
plot( steer(:,1), u);
hold on
plot( tempo, ones(size(tempo))*0.8e5)
plot( tempo, -ones(size(tempo))*0.8e5)
title('Rolagem')
xlabel('tempo (s)')
ylabel('\phi (rad)')