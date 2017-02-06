







%  p = plot(x,expAnalytical ); 
%   hold on
%  p1 = plot(x,expExplicitC025 ); 
%  hold on
%  p2 = plot(x,expExplicitC050 );
%   hold on
%  p3 = plot(x,expExplicitC075 );
%   hold on
%  p4 = plot(x,expExplicitC0999 );
% 
% xlabel('space[units]')
% ylabel('f(space. time)')
% title('Schemes compare results')
% legend('Analytical solution','Explicit Upwind Scheme', 'C =0.25' ,'C =0.5', 'C =0.75', 'C =0.999', 'C = 1.25', 'C =1.5', 'C =1.75', 'C =2' );
% %All C < 1
% legend('Analytical solution', 'C =0.25' ,'C =0.5', 'C =0.75', 'C =0.999');
% p(1).LineWidth = 2;
% p1(1).LineWidth = 2;
% p2(1).LineWidth = 2;
% p3(1).LineWidth = 2;
% p4(1).LineWidth = 2;
%% Different time

%  p = plot(x,expAnalytical ); 
%  hold on
%  p1 = plot(x,expExplicitTime10 );
%  hold on
%  xlabel('space[units]')
% ylabel('f(space. time)')
% title('Time compare results')
% legend('Time = 5', 'Time = 10');
% p(1).LineWidth = 2;
% p1(1).LineWidth = 2;

% figure
% waterfall(TimeAllTo10_2);

%% Different number of points
% 
%  p = plot(x,expExplicitTime10); 
%  hold on
% p1 = plot(x1,expExplicitPoints200);
% hold on
% p2 = plot(x3,expExplicitPoints400 );
% hold on
% xlabel('space[units]')
% ylabel('f(space. time)')
% title('Time compare results')
% legend('100 points', '200 points', '400 points');
% p(1).LineWidth = 2;
% p1(1).LineWidth = 2;
% p2(1).LineWidth = 2;