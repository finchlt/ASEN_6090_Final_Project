function [result] = readIgsTides(filename,reflla,dt)

if nargin == 1
	reflla=[];
elseif nargin == 2
	dt = [];
end

f=load(filename);

result.('t')=f(:,1)';

ne=length(result.t);

result.('env')=zeros(3,ne);
result.env(1,:)=f(:,2)';
result.env(2,:)=f(:,3)';
result.env(3,:)=f(:,4)';

if ~isempty(dt)
	tt = result.t(1):dt:result.t(end);
	e = interp1(result.t,result.env(1,:),tt,'linear');
	n = interp1(result.t,result.env(2,:),tt,'linear');
	u = interp1(result.t,result.env(3,:),tt,'linear');
	
	result.('t')=tt;
	
	ne=length(result.t);
	
	result.('env')=zeros(3,ne);
	result.env(1,:)=e;
	result.env(2,:)=n;
	result.env(3,:)=u;
end

if ~isempty(reflla)
	
	reflat = reflla(1);
	reflon = reflla(2);
	refalt = reflla(3);
	
	result.('xyz')=zeros(3,ne);
	
	for i=1:ne
		result.xyz(:,i) = enu2xyz(reflat, reflon, refalt, ...
			result.env(1,i),result.env(2,i),result.env(3,i));
	end
	
end

