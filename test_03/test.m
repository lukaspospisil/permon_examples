clear all

m = 15;

t = 0:2*pi/m:(2*pi-2*pi/m);

S = [-2,0; 2,0]';

P_points = zeros(2,length(t)); % the points of the first polytope
Q_points = zeros(2,length(t)); % the points of the second polytope

for i = 1:length(t)
    xP = [cos(t(i));sin(t(i))];
    xQ = [cos(pi - t(i));sin(pi - t(i))];
    
    P_points(:,i) = xP + S(:,1);
    Q_points(:,i) = xQ + S(:,2);
    
end

C = [P_points,-Q_points];

A = 2*C'*C;
b = zeros(2*m,1);

lb = zeros(2*m,1);
B = kron(eye(2),ones(1,m));
c = ones(2,1);

y = quadprog(A,-b,[],[],B,c,lb,[],[]);

point1 = P_points*y(1:m);
point2 = Q_points*y(m+1:2*m);

distance = 0.25*y'*A*y;

figure
hold on
fill([P_points(1,1:end) P_points(1,1)],[P_points(2,1:end) P_points(2,1)],[0.5 1 0.5]);
fill([Q_points(1,1:end) Q_points(1,1)],[Q_points(2,1:end) Q_points(2,1)],[0.5 0.5 1]);

plot([point1(1) point2(1)],[point1(2) point2(2)],'r');

my_legend = legend('polytope $\mathcal{P}$','polytope $\mathcal{Q}$');
set(my_legend,'Interpreter','latex')

axis equal
axis([-3.5 3.5 -1.5 2.0])
hold off

