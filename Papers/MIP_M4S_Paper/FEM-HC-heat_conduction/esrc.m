function y=esrc(x)
% 1MW, 90 fuel elements, Zr inner rod radius=0.003175m
% fuel meat radius=0.0174115m, fuel height = 0.381m
% 1e6/90/(pi*(0.0174115^2-0.003175^2)*.381)=3.1674e7 W/m3
% With ICONE document, we transform the heat flux value to the power
% density value, see source excel document : 2.84188e7 W/m3
y=3.1674e7; % W/m3
% y=3.1674e7*1.41966; % for comparison with Chance's thesis
end