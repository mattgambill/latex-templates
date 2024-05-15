function p = rpulse(t,tau)
% p = rpulse(t,tau)
% Creates a rectangular pulse of width tau from 0<t<tau
% t: the current time
% tau: the pulse width

    p = (t<=tau)&(t>=0);
end