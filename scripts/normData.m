function [Uj, pj, rhoj] = normData(test_name)
% input:
% folder:   path to folder of test case being analyzed

% output:
% Uj:   ideally expanded jet velocity under isentropic expansion
% pj:   ideally expanded jet pressure under isentropic expansion
% rhoj: ideally expanded jet density under isentropic expansion

% extract mach number
if contains(test_name, "M0p5"); Mj = 0.5; end
if contains(test_name, "M0p8"); Mj = 0.8; end
if contains(test_name, "M0p9"); Mj = 0.9; end

% assumptions
Tt = 1; gamma = 1.4; R = 1/gamma;
pj = 1/1.4;  % for all cases

% find Tj and cj
Tj = Tt * (1 + (gamma-1)/2 * Mj^2)^(-1);
cj = sqrt(Tj);

% use cj to find Uj
Uj = Mj * cj;

% use pj and Tj to find rhoj
rhoj = pj / (R * Tj);
end