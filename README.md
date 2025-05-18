# Calculating-Vth-change-of-ISFET-from-Id-change
The dataset and program are related to the work published as "Integration of AC Bias-Stabilized SiNRFETs with Microscale Ag/AgCl Pseudo-Reference Electrodes for Lab-on-chip Biosensing Application"
The example dataset was collected using ISFET sensors fabricated in our research group. The sensors were used to monitor the acidification of LB medium (supplemented with 1 wt% glucose) caused by the metabolic activity of E. coli.

sampling.txt records the drain current (ID) variations over time from multiple independent ISFET devices, corresponding to the gradual acidification of the culture medium.

IV.txt contains the ID–VG transfer curves of the ISFET sensors, used for threshold voltage (Vth) extraction.
Time.txt is defined the starting point of the monitoring. For example, Input 180 in time. txt. means the starting point of the monitoring is 180s.

All data are recorded using Hp4155A during experiment.
By processing this dataset with the provided script, a file named Vt1.xlsx is generated, which records the threshold voltage (Vth) as a function of time.
The Vth shift directly reflects pH changes in the culture medium, with a sensitivity of approximately 60 mV per pH unit.
