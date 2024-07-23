% -----------------------------------------------------------------------------
% Copyright (C) 2024 Marco Cattaneo
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of  MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% -----------------------------------------------------------------------------

function [t, up, lo] = ExperimentalPressureEnvelope

% Load experimentally recorded ultrasound pulse
data = load('UltrasoundPulse.dat');
fs = 1e10;
dt = 1/fs;

t = data(:,1);
P = data(:,2);

% Low-pass filter the pressure signal
fpass = 2e6;
P = lowpass(P,fpass,fs, 'ImpulseResponse','iir');

% Extract the upper and lower envelope
[up,lo] = envelope(P,5000,'peak');
env = 0.5*(up-lo);

% Normalisation envelopes
index_US_start = find(env>0.00005, 1, 'first');
t = dt*(0:length(t(index_US_start:end))-1);
scalefactor = 1/max(up);
up = scalefactor*up(index_US_start:end);
lo = -scalefactor*lo(index_US_start:end);






