function [ distance ] = dist( v1,v2 )

if ( (size(v2,2) == size(v1,1)) && (size(v2,1) == size(v1,2) ) )
	v2 = v2';
end
v = v2-v1;
%distance = sqrt( (v2(1)-v1(1))^2 + (v2(2)-v1(2))^2 + (v2(3)-v1(3))^2 );
distance = sqrt(v(1)^2 + v(2)^2 + v(3)^2);
