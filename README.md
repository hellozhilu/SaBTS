# SaBTS
Source codes for the article: "Stagnation-aware breakout tabu search for the minimum conductance graph partitioning problem"

Zhi Lu, Jin-Kao Hao and Yi Zhou. "Stagnation-aware breakout tabu search for the Minimum Conductance Graph Partitioning Problem"

The source codes of our proposed Greedy+SaBTS and Metis+SaBTS algorithm, along with the reimplemented StS-AMA algorithm are available in here. 
Please make sure that the above paper is cited if you use the code in your research. 
The software is distributed for academic purposes only. 
If you wish to use this software for commercial applications, please obtain permission from Zhi Lu (zhilusix@gmail.com) or Jin-Kao Hao (jin-kao.hao@univ-angers.fr).


---------------------------------------------------------------------------------------------------------------------------------------------------------------------


How to run the programs

Greedy+SaBTS:
./${ExecuteFile}" ${InstanceFile} ${DataSet} ${RunTime} ${NumberRepeats} ${alfa} ${D} ${T} ${L0} ${P0}

Metis+SaBTS:
./${ExecuteFile}" ${InstanceFile} ${MetisFile} ${DataSet} ${RunTime} ${NumberRepeats} ${alfa} ${D} ${T} ${L0} ${P0}

StS-AMA:
./${ExecuteFile}" ${InstanceFile} ${DataSet} ${RunTime} ${NumberRepeats}


// For example to run Greedy+SaBTS algorithm
// Given an instance, e.g., karate.graph
InstanceFile="karate.graph"
DataSet="test"
RunTime=3600
NumberRepeats=20
alfa=100
D=6000
T=1000
L0=0.4
P0=0.8
