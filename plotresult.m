set(0, 'DefaultLineLineWidth', 2);
subplot(2,2,1)
tempo = steer(:,1);
plot( tempo, x(1,:)'*180/pi);
hold on
if ~isempty(ub)
    plot( tempo, ones(size(tempo))*ub*180/pi)
    plot( tempo, ones(size(tempo))*lb*180/pi)
end
title('Deslizamento lateral')
xlabel('tempo (s)')
ylabel('\beta (graus)')

subplot(2,2,2)
plot( steer(:,1), x(2,:)'*180/pi);
hold on
plot( steer(:,1), desired_yaw_rate*180/pi);
%plot( tempo, ones(size(tempo))*pi/6)
%plot( tempo, -ones(size(tempo))*pi/6)
title('Taxa de guinada')
xlabel('tempo (s)')
ylabel('d\psi/dt (graus/s)')

subplot(2,2,3)
plot( steer(:,1), x(4,:)'*180/pi);
hold on
%plot( tempo, ones(size(tempo))*pi/18)
%plot( tempo, -ones(size(tempo))*pi/18)
title('Rolagem')
xlabel('tempo (s)')
ylabel('\phi (graus)')

subplot(2,2,4)
plot( steer(:,1), u);
hold on
if ~isempty( cmd_ub )
    plot( tempo, ones(size(tempo))*cmd_ub)
    plot( tempo, ones(size(tempo))*cmd_lb)
end
title('Comando')
xlabel('tempo (s)')
ylabel('M_u (Nm)')
