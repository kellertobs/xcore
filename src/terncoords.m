% TERNCOORDS calculate rectangular coordinates of fractions on a ternary plot
%   [X, Y] = TERNCOORDS(FA, FB) returns the rectangular X and Y coordinates
%   for the point with a fraction defined by FA and FB.  It is assumed that
%   FA and FB are sensible fractions.
%
%   [X, Y] = TERNCOORDS(FA, FB, FC) returns the same.  FC is assumed to be
%   the remainder when subtracting FA and FB from 1.

% Author: Carl Sandrock 20050211

% Modified: Tobias Keller 20250214

function [x, y] = terncoords(fA, fB, fC)

if nargin < 3 && any(fA+fB>1)
    error('Invalid data input: normalise end-member fractions to unit sum or provide all three end-members!');
elseif nargin < 3
    fC = 1 - (fA + fB);
else
    sumf = (fA + fB + fC);
    fA = fA./sumf;
    fB = fB./sumf;
    fC = fC./sumf;
end

y = fB*sin(deg2rad(60));
x = fA + y*cot(deg2rad(60));
