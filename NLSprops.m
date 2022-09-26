classdef NLSprops
   properties (Constant)
      Nw = 2 % window length [T1]
      Ns = 0.05 % sliding length [T1]
      freqdev = 0.2 % frequency variation allowed for, fraction
      ampdev = 2 % amplitude variation allowed for, fraction
      dampthresh = 0.9 % damping ratio elimination criterion, [-]
      plotflag = 1 % plotting option
      drawinstant = 1 % plotting option
   end
end
