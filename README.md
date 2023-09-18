# eegstudy
TFCE.mat and TFCE.mat is the mat.file store the MSE and Relative power spectrum results separately.
It contains: AAopenclosed_lastordered; CAclosed_lastordered; CAopen_lastordered; CCopenclosed_lastordered and DiffAC_lastordered, represents for:
two sample t test with contrast between eyes-open and eyes-closed condition in ASC group; contrasts between control and ASC group at the eyes-closed contion and at the eyes-open condition;
and contrasts between eyes-open and eye-closed condition in control group.

In each variable, it contains p value, raw t statistics and TFCE statistics for each channel and time scales, resulting in 32X20 matrix.

tfce_script is for data input and heatmap output for each condition.
