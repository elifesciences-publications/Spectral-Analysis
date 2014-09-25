

[master.Data, master.Time,master.Time_reduced, master.statMat] =...
    deal(dataStruct,timeAxisStruct,time_reduced,statMat);

[master.Wxy, master.Wxy_avg, master.meanPhaseVec, master.maxPhaseVec] =...
    deal(Wxy, Wxy_avg, meanPhaseVec, maxPhaseVec);

[master.meanPowVec, master.meanFreqVec, master.maxFreqVec] = ...
    deal(tvpower, mfvec, pfvec);




