function [h, m, s] = sec2hms(sec)
h = 0;
m = 0;

if sec > 3600
	h = floor(sec/3600);
end

if sec>60
	m = floor((sec - 3600*h)/60);
end

s = sec - m*60 - h*3600;

end