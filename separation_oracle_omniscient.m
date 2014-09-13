function [c,gamma] = separation_oracle_omniscient(C, d, a)
% SEPARATION_ORACLE_OMNISCIENT Helper function for the ellipsoid method.
%   [c,gamma] = SEPARATION_ORACLE_OMNISCIENT(C,d,a) will check any
%   inequality in C,d is violated by a and return the separation hyperplane
%   c, otherwise it will return a as cutting plane.
%
%   See also ellipsoid_method

%% Separation oracle
% This separation oracle is called with all avaiable information
% (Constraint matrix C, Values of inequalities d, current center point a).
% If all constraints are satisfied, point a is in the polytope given
% by C and it will be returned. Otherwise a violated constraint is returned
% as cutting plane. Gamma is the value of the violated constraint, which is
% need for calculating the deep cut.
%
% * Author: Florian Stallmann <mail@florian-stallmann.de>
% * Licence: CC0 1.0 Universal (Public Domain Dedication)
% * Version: 2014-09-08

% Check if all constraints are satisfied
if all(C*a <= d)
    c = a; % a is in polytope, so return a
    gamma = 0; % gamma not needed
else
    % find any violated constraint
    j = find(C*a > d, 1);
    % return violated inequality as separating hyperplane
    c = C(j,:)';
    % value of inequality (needed for deep cut)
    gamma = d(j);
end