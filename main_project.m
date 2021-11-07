clear all; close all;

zz = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];

h1 = ez;
h2 = ey;
h3 = ey;
h4 = ey;

L1 = 100;
L2 = 135;
L3 = 160;
L4 = 50;
L5 = 75;

p01 = L1*ez;
p12 = zz;
p23 = L2*ex;
p34 = L3*ex;
p4T = L4*ex-L5*ez;

dobot.P=[p01 p12 p23 p34 p4T];
dobot.H=[h1 h2 h3 h4];
dobot.joint_type=[0 0 0 0];
dobot.q=[0;0;0;0];

radius=.01;
[dobot_rbt,colLink]=defineRobot(dobot,radius);

qrs = zeros(1, 4);
qrs(1:3) = (pi()/2 + pi()/2).*rand(3,1) - pi()/2;
qrs(4) = -1*(qrs(2)+qrs(3));
dobot.q = qrs;
dobot=fwddiffkiniter(dobot);
dobot=invkin(dobot);
for i=1:4
    if abs(qrs-dobot.q(:,i)) < 1e-4
        disp('Random q');
        disp(qrs);
        disp('Calculated q')
        disp(dobot.q(:,i));
    end   
end

function robot = invkin(robot)

zz = [0;0;0];
ex = [1;0;0];
ey = [0;1;0];
ez = [0;0;1];
h1=robot.H(:,1);
h2=robot.H(:,2);
h3=robot.H(:,3);
h4=robot.H(:,4);

qsol1 = [0;0;0;0];
qsol2 = [0;0;0;0];
qsol3 = [0;0;0;0];
qsol4 = [0;0;0;0];

T = robot.T;
R0T = T(1:3,1:3);
p0T = T(1:3, 4);

p01 = robot.P(:,1);
p12 = robot.P(:,2);
p23 = robot.P(:,3);
p34 = robot.P(:,4);
p4T = robot.P(:,5);

p14_0 = p0T - p01 - R0T*p4T;

%% solve for q3
q3 = subprob3(h3, -p34, p23, norm(p14_0));

qsol1(3) = q3(1);
qsol2(3) = q3(1);
qsol3(3) = q3(2);
qsol4(3) = q3(2);


%% solve for q1 and q2
[q1_1,q2_1]=subprob2(-h1,h2,p14_0,p23+rot(h3, q3(1))*p34);
[q1_2,q2_2]=subprob2(-h1,h2,p14_0,p23+rot(h3, q3(2))*p34);
qsol1(1:2) = [q1_1(1); q2_1(1)];
qsol2(1:2) = [q1_1(2); q2_1(2)];
qsol3(1:2) = [q1_2(1); q2_2(1)];
qsol4(1:2) = [q1_2(2); q2_2(2)];

%% solve for q4 

qsol1(4) = -1*(qsol1(2)+qsol1(3));
qsol2(4) = -1*(qsol2(2)+qsol2(3));
qsol3(4) = -1*(qsol3(2)+qsol3(3));
qsol4(4) = -1*(qsol4(2)+qsol4(3));

qsol = [qsol1 qsol2 qsol3 qsol4];

robot.q = qsol;

end

function robot=fwddiffkiniter(robot)

q=robot.q;
n=length(robot.q);

T=eye(4,4);
joint_type=robot.joint_type;
if ~exist('joint_type');robot.joint_type=zeros(1,n);end

for i=1:n
    h=robot.H(1:3,i);
    if robot.joint_type(i)==0
        R=expm(hat(h)*q(i));p=robot.P(1:3,i);
    else
        R=eye(3,3);p=robot.P(1:3,i)+q(i)*h;
    end
    T=T*[R p;zeros(1,3) 1];
end    
    
robot.T=T*[eye(3,3) robot.P(1:3,n+1);zeros(1,3) 1];

end
%% Subproblems

%
% q=subprob3(k,p1,p2,d)
%
% solve for theta from
%
% norm(p2-exp(k x q) p1) = d
%
% input: k,p1,p2 as R^3 column vectors, delta: scalar
% output: q (2x1 vector, 2 solutions)
%

function q=subprob3(k,p1,p2,d)

pp1=p1-k'*p1*k;
pp2=p2-k'*p2*k;
dpsq=d^2-(k'*(p1-p2))^2;

if dpsq<0
    theta=[NaN;NaN];
    return;
end

if dpsq==0
    theta=subprob1(k,pp1/norm(pp1),pp2/norm(pp2));
    return;
end
 
bb=(norm(pp1)^2+norm(pp2)^2-dpsq)/(2*norm(pp1)*norm(pp2));
if abs(bb)>1
    theta=[NaN;NaN];
    return
end

phi=acos(bb);

q0=subprob1(k,pp1/norm(pp1),pp2/norm(pp2));
q=zeros(2,1);

q(1)=q0+phi;
q(2)=q0-phi;
end

%
% [q1,q2]=subprob2(k1,k2,p1,p2)
%
% solve for theta1 and theta2 from
%
% exp(k1 x q1) p1 = exp(k2 x q2) p2 
%
% input: k1,k2,p1,p2 as R^3 column vectors
%
% output: q1 and q2 as 2x1 columns corresponding to the two solutions
%

function [q1,q2]=subprob2(k1,k2,p1,p2)

p2=p2/norm(p2)*norm(p1);
k12=k1'*k2;
pk1=p1'*k1;
pk2=p2'*k2;

% check if solution exists

if abs(k12^2-1)<eps;theta1=[];theta2=[];
    q1=[NaN;NaN];q2=[NaN;NaN];
    disp('no solution (k1 and k2 are collinear)');
    return;
end

a=[1 -k12; -k12 1]*[pk1;pk2]/(1-k12^2);

% 
% check if solution exists
%
cond=(norm(p1)^2-norm(a)^2-2*a(1)*a(2)*k12);

% special case: 1 solution
if abs(cond)<eps;
  v=[k1 k2]*a;
  q1a=subprob1(k1,p1,v);
  q2a=subprob1(k2,p2,v);
  q1=[q1a;q1a];
  q2=[q2a;q2a];
end

% special case: no solution
if cond<0
    q1=[NaN NaN];q2=[NaN NaN];
    disp('no solution (two cones do not intersect)');
    return;
end

gamma=sqrt(cond)/norm(cross(k1,k2));

% general case: 2 solutions

q1=zeros(2,1);
q2=zeros(2,1);

v1=[k1 k2 cross(k1,k2)]*[a;gamma];
v2=[k1 k2 cross(k1,k2)]*[a;-gamma];
q1(1)=subprob1(k1,p1,v1);
q1(2)=subprob1(k1,p1,v2);

q2(1)=subprob1(k2,p2,v1);
q2(2)=subprob1(k2,p2,v2);

end

%
% q=subprob1(k,p1,p2)
%
% solve for q from
%
% exp(k x q) p1 = p2
%
% input: k,p1,p2 as R^3 column vectors
% output: q (scalar)
%

function q=subprob1(k,p1,p2)

p2=p2/norm(p2)*norm(p1);

if norm(p1-p2)<sqrt(eps);q=0;return;end
  
k=k/norm(k);
pp1=p1-(p1'*k)*k;
pp2=p2-(p2'*k)*k;

epp1=pp1/norm(pp1);
epp2=pp2/norm(pp2);

q=subprob0(k,epp1,epp2);
%q=atan2(k'*(cross(epp1,epp2)),epp1'*epp2);
end

%
% q=subprob0(k,p1,p2)
%
% solve for q subtended between p1 and p2
%    k determines the sign of q
%
% input: k,p1,p2 as R^3 column vectors
% output: q (scalar)
%

function q=subprob0(k,p1,p2)

if ((k'*p1)>sqrt(eps)|(k'*p2)>sqrt(eps))
  error('k must be perpendicular to p and q');
end

p1=p1/norm(p1);
p2=p2/norm(p2);

q=2*atan2(norm(p1-p2),norm(p1+p2));

if k'*(cross(p1,p2))<0
  q=-q;
end
end

function [robot,colLink]=defineRobot(robdef,rad)

% product of exponential description
H=robdef.H;
P=robdef.P;
type=robdef.joint_type;
n=numel(type);

% homogeneous transforms T0i stored in T0i{k}, k=1,..,n+1
T=eye(4,4);
for i=1:n
    q(i)=0;
    if type(i)==0
        Ti{i}=[rot(H(:,i),q(i)) P(:,i);[0 0 0 1]];
        T=T*Ti{i};
    else
        Ti{i}=[eye(3,3),P(:,i)+q(i)*H(:,i);[0 0 0 1]];
        T=T*Ti{i};
    end
    T0i{i}=T;
end
Ti{n+1}=[eye(3,3) P(:,n+1);[0 0 0 1]];
T0i{n+1}=T*Ti{n+1};

% define MATLAB rigidbody tree

robot=rigidBodyTree('DataFormat','column');

for i=1:n
    ii=num2str(i);
    eval(['body',ii,' = rigidBody(''','body',ii,''');']);
    if type(i)==0
        eval(['jnt',ii,' = rigidBodyJoint(''','jnt',...
            ii,'''',',''revolute''',');']);
    else
        eval(['jnt',ii,' = rigidBodyJoint(''','jnt',...
            ii,'''',',''prismatic''',');']);
    end
    eval(['jnt',ii,'.JointAxis=H(:,',ii,');']);
end
ii=num2str(n+1);
eval(['body',ii,' = rigidBody(''','body',ii,''');']);
eval(['jnt',ii,' = rigidBodyJoint(''','jnt',...
    ii,'''',',''fixed''',');']);

for i=1:n+1
    ii=num2str(i);
    eval(['setFixedTransform(jnt',ii,',Ti{',ii,'});']);
    eval(['body',ii,'.Joint = jnt',ii,';']);
end

addBody(robot,body1,'base')
for i=2:n+1
    ii=num2str(i);
    iim1=num2str(i-1);
    eval(['addBody(robot,body',ii,',','''body',iim1,''');']);
end

colBodyRadius=rad;

for i=1:n+1
    if norm(P(:,i))>0
        if (i<n+1) && (type(i)>0)
            colLink{i} = collisionBox(colBodyRadius,colBodyRadius,norm(P(:,i)));
        else
            colLink{i} = collisionCylinder(colBodyRadius,norm(P(:,i)));
        end
        kvec=cross([0;0;1],P(:,i));
        if norm(kvec)<sqrt(eps)
            colLink{i}.Pose=trvec2tform(P(:,i)'/2);
        else
            th=subprob0(kvec,[0;0;1],P(:,i));
            colLink{i}.Pose=[rot(kvec,th) P(:,i)/2;[0 0 0 1]];
        end
    else
        colLink{i} = collisionCylinder(colBodyRadius,norm(P(:,i)));
    end
    if i==1
        addCollision(robot.Base,colLink{i});    
    else
        addCollision(robot.Bodies{i-1},colLink{i});    
    end    
end

end