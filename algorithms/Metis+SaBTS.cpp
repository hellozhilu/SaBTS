//============================================================================
// Name        : Metis+SaBTS.cpp
// Author      : Zhi LU (zhilusix@gmail.com)
// Version     : Feb. 2018
// Copyright   : LERIA, Faculté d'Sciences, Université d'Angers, France.
// Description : Stagnation-awared Breakout Tabu Search
//               for the Minimum Conductance Graph Partitioning Problem.
//               WITH *Metis initialization*
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define MAX_INT 10000000000
#define EPSILON 0.000001
#define MAX(X,Y)  ((X) >= (Y) ? (X) : (Y))
#define MIN(X,Y)  ((X) <= (Y) ? (X) : (Y))
#define ABS(X)  (((X) < 0) ? -(X) : (X))

using namespace std;


/* global variable */
//graph
long long int nbNodes;
long long int nbEdges;
long long int *adjacLength;
long long int **adjacTable;
double graphDensity;

long long int volumeS, volumeVS, cG;  		//cG: cutting edges number
long long int sizeS, sizeVS, sizeCut;
long long int localBestSol, globalBestSol, metisBestSol;
int *setFlag, *globalsetFlag,  *metisFlag;  //labeled nodes belong to which sets, 0-S, 1-V\S
long long int *degvtoS, *degvtoVS;
long long int *setVS;
long long int *setCut, *posCut;      	    //setCut[]: stored nodes related to the cut edge
int *weiNodes;                       		//vertex weighting array for vertex weighting scheme
int *tabuTenure;
int w;                                      //counter of non-improving local optima
int L;                                      //jump magnitude
//parameters (5)
int alfa;                                   //tabu tenure management factor = 100
int D;                                      //depth of tabu search = 6000
int T;                                      //stagnation threshold = 1000
double L0;                                  //initial jump magnitude = 0.4 * |V|
double P0;                       	        //minimum probability = 0.8

//record time
int runTime;
int nbRepeat;
struct timespec startTime = {0, 0};
struct timespec endTime = {0, 0};
double globalTime;
long long int seed;


/* function prototype */
void read_graph(char*);
void read_metis(char*);
void allocate_memory();
void check_sol(int*, long long int&);
void check_result();
double diff_time(timespec);
double rand_double();
void set_node_weight(long long int&);
void initialize(int*, long long int&);
void initialize_metis_input();
void cut_set_update(int*, long long int&);
void cut_neighboring_tabu_search();
void quick_sort(int*, int*, int, int);
void local_optimal_exploring_ADD(int&);
void local_optimal_exploring_SWAP(int&);
void local_optimal_escaping(int&);
void free_memory();
void SaBTS();


/* read instances */
void read_graph(char *filename)
{
	string line;
	ifstream ff(filename);
	if (ff.fail())
	{
		fprintf(stderr, "Error in read_graph()-1\n");
		exit(-1);
	}
	getline(ff, line);
	while ((line.length()==0) || (line[0]=='%'))
		getline(ff, line);

	stringstream ss(line);

	int nbType = 0;
	ss >> nbNodes >> nbEdges >> nbType;

	//graph density
	long long int completeEdges;
	completeEdges = nbNodes*(nbNodes-1)/2;
	graphDensity = (double)nbEdges/(double)completeEdges;

	adjacTable = new long long int*[nbNodes];
	adjacLength = new long long int[nbNodes];
	long long int *tempList = new long long int[nbNodes];
	long long int len = 0, lno = 0, vtx = 0;

	while (getline(ff, line))
	{
		if (ff.fail())
		{
			fprintf(stderr, "Error in read_graph()-2\n");
			exit(-1);
		}

		len = 0;
		stringstream ss(line);

		while (ss >> vtx)
			tempList[len++] = vtx-1;
		adjacLength[lno] = len;

		if (len == 0)
			adjacTable[lno] = NULL;
		else
		{
			adjacTable[lno] = new long long int[len];
			memcpy(adjacTable[lno], tempList, sizeof(long long int)*len);
		}

		lno++;
	}

	//allocate memory
	allocate_memory();

	//free memory
	delete[] tempList;
	tempList = NULL;

	cout << "finish reading" << endl;
	cout << "nbNodes=" << nbNodes << ", nbEdges=" << nbEdges << ", graphDensity=" << graphDensity << endl << endl;

	ff.close();
}


void read_metis(char *filename)
{
	string line;
	ifstream ff(filename);
	if (ff.fail())
	{
		fprintf(stderr, "Error in read_metis()-1\n");
		exit(-1);
	}

	metisFlag = new int[nbNodes];
	int i = 0;

	while (!ff.eof())
	{
		getline(ff, line);
		stringstream ss(line);

		if (line[0] == '0')
			metisFlag[i] = 0;
		else if (line[0] == '1')
			metisFlag[i] = 1;

		i++;
	}

	ff.close();
}


void allocate_memory()
{
	setFlag = new int[nbNodes];
	globalsetFlag = new int[nbNodes];
	degvtoS = new long long int[nbNodes];
	degvtoVS = new long long int[nbNodes];
	setVS = new long long int[nbNodes];
	setCut = new long long int[nbNodes];
	posCut = new long long int[nbNodes];
	weiNodes = new int[nbNodes];
	tabuTenure = new int[nbNodes];
}


/* verify intermediate result */
void check_sol(int *set, long long int &sol)
{
	//1.check volumeS, volumeVS, cG
	long long int volumeS_check = 0, volumeVS_check = 0, cG_check = 0;
	for (long long int i=0; i<nbNodes; i++)
	{
		if (set[i] == 0)
		{
			volumeS_check += adjacLength[i];
			for (long long int j=0; j<adjacLength[i]; j++)
				if (set[adjacTable[i][j]] == 1)
					cG_check += 1;
		}
		else
			volumeVS_check += adjacLength[i];
	}
	if (volumeS != volumeS_check)
	{
		printf("volumeS=%lld, volumeS_check=%lld\n", volumeS, volumeS_check);
		fprintf(stderr, "Error volumeS in check_sol()\n");
		exit(-1);
	}
	if (volumeVS != volumeVS_check)
	{
		printf("volumeVS=%lld, volumeVS_check=%lld\n", volumeS, volumeS_check);
		fprintf(stderr, "Error volumeVS in check_sol()\n");
		exit(-1);
	}
	if (cG != cG_check)
	{
		printf("cG=%lld, cG_check=%lld\n", cG, cG_check);
		fprintf(stderr, "Error cG in check_sol()\n");
		exit(-1);
	}

	//2.check degvtoS[], degvtoVS[]
	long long int *degvtoS_check;
	long long int *degvtoVS_check;
	degvtoS_check = new long long int[nbNodes];
	degvtoVS_check = new long long int[nbNodes];
	for (long long int i=0; i<nbNodes; i++)
	{
		degvtoS_check[i] = 0;
		degvtoVS_check[i] = 0;
		for (long long int j=0; j<adjacLength[i]; j++)
			if (set[adjacTable[i][j]] == 0)
				degvtoS_check[i]++;
			else
				degvtoVS_check[i]++;
	}
	for (long long int i=0; i<nbNodes; i++)
	{
		if (degvtoS[i] != degvtoS_check[i])
		{
			printf("Error no. %lld\n", i);
			fprintf(stderr, "Error degvtoS[] in check_sol()\n");
			exit(-1);
		}
		if (degvtoVS[i] != degvtoVS_check[i])
		{
			printf("Error no. %lld\n", i);
			fprintf(stderr, "Error degvtoVS[] in check_sol()\n");
			exit(-1);
		}
	}

	//3.check objective function
	long long int sol_check;
	sol_check = (double)cG_check/(double)MIN(volumeS_check, volumeVS_check)*MAX_INT;
	if (sol != sol_check)
	{
		printf("sol=%.10f, sol_check=%.10f\n", (double)sol/MAX_INT, (double)sol_check/MAX_INT);
		fprintf(stderr, "Error sol in check_sol()\n");
		exit(-1);
	}

	//4.check sizeVS, setVS[]
	//sort
	long long int t;
	for (long long int i=0; i<sizeVS; i++)
	{
		for (long long int j=i+1; j<sizeVS; j++)
		{
			if (setVS[i] > setVS[j])
			{
				t = setVS[i];
				setVS[i] = setVS[j];
				setVS[j] = t;
			}
		}
	}
	//check
	long long int sizeVS_check = 0;
	long long int *setVS_check;
	setVS_check = new long long int[nbNodes];
	for (long long int i=0; i<nbNodes; i++)
		if (setFlag[i] == 1)
			setVS_check[sizeVS_check++] = i;
	if (sizeVS != sizeVS_check)
	{
		printf("sizeVS=%lld, sizeVS_check=%lld\n", sizeVS, sizeVS_check);
		fprintf(stderr, "Error sizeVS in check_sol()\n");
		exit(-1);
	}
	for (long long int i=0; i<sizeVS_check; i++)
	{
		if (setVS[i] != setVS_check[i])
		{
			printf("Error no. %lld\n", i);
			fprintf(stderr, "Error setVS[] in check_sol()\n");
			exit(-1);
		}
	}

	//5.check setCut[], posCut[] and sizeCut
	//check1
	for (long long int i=0; i<sizeCut; i++)
	{
		if (posCut[setCut[i]] != i)
		{
			printf("Error no. %lld\n", i);
			fprintf(stderr, "Error cutSet[] and posCut[] in check_sol()\n");
			exit(-1);
		}
	}
	//check2
	//sort
	long long int k;
	for (long long int i=0; i<sizeCut; i++)
	{
		for (long long int j=i+1; j<sizeCut; j++)
		{
			if (setCut[i] > setCut[j])
			{
				k = setCut[i];
				setCut[i] = setCut[j];
				setCut[j] = k;
			}
		}
	}
	for (long long int i=0; i<nbNodes; i++)
		posCut[i] = -1;
	for (long long int i=0; i<sizeCut; i++)
		posCut[setCut[i]] = i;
	//check
	long long int sizeCut_check = 0;
	long long int *setCut_check;
	long long int *posCut_check;
	setCut_check = new long long int[nbNodes];
	posCut_check = new long long int[nbNodes];
	for (long int i=0; i<nbNodes; i++)
		posCut_check[i] = -1;
	for (long long int i=0; i<nbNodes; i++)
	{
		for (long long int j=0; j<adjacLength[i]; j++)
		{
			if (set[i] != set[adjacTable[i][j]])
			{
				posCut_check[i] = sizeCut_check;
				setCut_check[sizeCut_check++] = i;
				break;
			}
		}
	}
	for (long long int i=0; i<nbNodes; i++)
	{
		if (posCut[i] != posCut_check[i])
		{
			printf("Error no. %lld\n", i);
			fprintf(stderr, "Error posCut[] in check_sol()\n");
			exit(-1);
		}
	}
	for (long long int i=0; i<sizeCut; i++)
	{
		if (setCut[i] != setCut_check[i])
		{
			printf("Error no. %lld\n", i);
			fprintf(stderr, "Error cutSet[] in check_sol()\n");
			exit(-1);
		}
	}

	//free memory
	delete[] degvtoS_check;
	delete[] degvtoVS_check;
	delete[] setVS_check;
	delete[] setCut_check;
	delete[] posCut_check;
	degvtoS_check = NULL;
	degvtoVS_check = NULL;
	setVS_check = NULL;
	setCut_check = NULL;
	posCut_check = NULL;
}


/* verify final result */
void check_result()
{
	long long int min_sol;

	volumeS = 0, volumeVS = 0, sizeS = 0, sizeVS = 0, cG = 0;
	for (long long int i=0; i<nbNodes; i++)
	{
		if (globalsetFlag[i] == 0)
		{
			sizeS++;
			volumeS += adjacLength[i];
			for (long long int j=0; j<adjacLength[i]; j++)
				if (globalsetFlag[adjacTable[i][j]] == 1)
					cG += 1;
		}
		else
		{
			sizeVS++;
			volumeVS += adjacLength[i];
		}
	}
	min_sol = (double)cG/(double)MIN(volumeS, volumeVS)*MAX_INT;

	if (ABS(min_sol-globalBestSol) > EPSILON)
	{
		printf("Find a error solution!\n");
		printf("min_sol=%.10f", (double)min_sol/(double)MAX_INT);
		printf(", while globalBestSol=%.10f\n", (double)globalBestSol/(double)MAX_INT);
		exit(0);
	}
}


/* record time */
double diff_time(timespec start)
{
	timespec temp;
	double spentTime;

	clock_gettime(CLOCK_REALTIME, &endTime);

	if ((endTime.tv_nsec-start.tv_nsec) < 0)
	{
		temp.tv_sec = endTime.tv_sec - start.tv_sec - 1;
		temp.tv_nsec = 1000000000 + endTime.tv_nsec - start.tv_nsec;
	}
	else
	{
		temp.tv_sec = endTime.tv_sec - start.tv_sec;
		temp.tv_nsec = endTime.tv_nsec - start.tv_nsec;
	}

//	spentTime = (double)(temp.tv_sec*1000) + (double)(temp.tv_nsec/1000000);                   //for ms
	spentTime = (double)temp.tv_sec + (double)(temp.tv_nsec/1000)/(double)CLOCKS_PER_SEC;      //for s

	return spentTime;
}


/* generate random decimal from (0,1) */
double rand_double()
{
    double r;

    r = (double)rand()/(double)RAND_MAX;
//	r = (rand()%10000)*0.0001;

	return r;
}


/* node weighting scheme */
void set_node_weight(long long int &v)
{
	for (long long int i=0; i<nbNodes; i++)
		weiNodes[i]++;

	weiNodes[v] = 0;
}


/* initialization */
void initialize(int *set, long long int &sol)
{
	//setV_S[], volumeS, volumeVS, cG, sizeS, sizeVS
	volumeS = 0, volumeVS = 0, cG = 0;
	sizeS = 0, sizeVS = 0;
	for (long long int i=0; i<nbNodes; i++)
	{
		if (set[i] == 0)
		{
			volumeS += adjacLength[i];
			sizeS++;
			for (long long int j=0; j<adjacLength[i]; j++)
				if (set[adjacTable[i][j]] == 1)
					cG += 1;
		}
		else
		{
			volumeVS += adjacLength[i];
			setVS[sizeVS++] = i;
		}
	}

	//ensure the solution is legal
	if ((volumeS==0) || (volumeVS==0))
	{
		printf("volumeS=%lld, volumeVS=%lld\n", volumeS, volumeVS);
		fprintf(stderr, "Error sol in initialize()1\n");
		exit(-1);
	}

	//degvtoS[], degvtoVS[]
	for (long long int i=0; i<nbNodes; i++)
	{
		degvtoS[i] = 0;
		degvtoVS[i] = 0;
		for (long long int j=0; j<adjacLength[i]; j++)
		{
			if (set[adjacTable[i][j]] == 0)
				degvtoS[i]++;
			else
				degvtoVS[i]++;
		}
	}

	//cutSet[], posCut[], sizeCut
	sizeCut = 0;
	for (long long int i=0; i<nbNodes; i++)
	{
		posCut[i] = -1;
		for (long long int j=0; j<adjacLength[i]; j++)
		{
			if (set[i] != set[adjacTable[i][j]])
			{
				setCut[sizeCut] = i;
				posCut[i] = sizeCut++;
				break;
			}
		}
	}

	//localBestSol
	long long int check_sol;
	check_sol = (double) cG / (double) MIN(volumeS, volumeVS) * MAX_INT;
	if (sol != check_sol)
	{
		printf("sol=%.10lf, check_sol=%.10lf\n", (double) sol / (double) MAX_INT, (double) check_sol / (double) MAX_INT);
		fprintf(stderr, "Error sol in initialize()2\n");
		exit(-1);
	}
}


void initialize_metis_input()
{
	for (long long int i = 0; i < nbNodes; i++)
		setFlag[i] = metisFlag[i];

	// initial variables and structures
	//setV_S[], volumeS, volumeVS, cG, sizeS, sizeVS
	volumeS = 0, volumeVS = 0, cG = 0;
	sizeS = 0, sizeVS = 0;
	for (long long int i = 0; i < nbNodes; i++)
	{
		if (setFlag[i] == 0)
		{
			volumeS += adjacLength[i];
			sizeS++;
			for (long long int j = 0; j < adjacLength[i]; j++)
				if (setFlag[adjacTable[i][j]] == 1)
					cG += 1;
		}
		else
		{
			volumeVS += adjacLength[i];
			setVS[sizeVS++] = i;
		}
	}

	//ensure the solution is legal
	if ((volumeS == 0) || (volumeVS == 0))
	{
		printf("volumeS=%lld, volumeVS=%lld\n", volumeS, volumeVS);
		fprintf(stderr, "Error sol in initialize()1\n");
		exit(-1);
	}

	//degvtoS[], degvtoVS[]
	for (long long int i = 0; i < nbNodes; i++)
	{
		degvtoS[i] = 0;
		degvtoVS[i] = 0;
		for (long long int j = 0; j < adjacLength[i]; j++)
		{
			if (setFlag[adjacTable[i][j]] == 0)
				degvtoS[i]++;
			else
				degvtoVS[i]++;
		}
	}

	//cutSet[], posCut[], sizeCut
	sizeCut = 0;
	for (long long int i = 0; i < nbNodes; i++)
	{
		posCut[i] = -1;
		for (long long int j = 0; j < adjacLength[i]; j++)
		{
			if (setFlag[i] != setFlag[adjacTable[i][j]])
			{
				setCut[sizeCut] = i;
				posCut[i] = sizeCut++;
				break;
			}
		}
	}

	localBestSol = (double) cG / (double) MIN(volumeS, volumeVS) * MAX_INT;
	metisBestSol = localBestSol;

	if (localBestSol < globalBestSol)
	{
		globalBestSol = localBestSol;
		globalTime = diff_time(startTime);
		for (long long int i = 0; i < nbNodes; i++)
			globalsetFlag[i] = setFlag[i];
		w=0;
	}
	else
		w++;
}


/* cut set updating */
void cut_set_update(int *set, long long int &move_node)
{
	long long int nb, k;

	nb = 0;
	//from V\S(1) to S(0)
	if(set[move_node] == 1)
	{
		//not in the cut set
		if (degvtoS[move_node] == 0)
		{
			//add itself
			setCut[sizeCut] = move_node;
			posCut[move_node] = sizeCut;
			sizeCut++;

			while (nb < adjacLength[move_node])
			{
				k = adjacTable[move_node][nb];
				if (degvtoVS[k] == 0)
				{
					fprintf(stderr, "Error in cut_set_update()\n");
					exit(-1);
				}
				if (degvtoS[k] == 0)
				{
					setCut[sizeCut] = k;
					posCut[k] = sizeCut;
					sizeCut++;
				}
				nb++;
			}
		}
		//in the cut set
		else
		{
			//consider itself
			if (degvtoVS[move_node] == 0)
			{
				setCut[posCut[move_node]] = setCut[sizeCut-1];
				posCut[setCut[sizeCut-1]] = posCut[move_node];
				posCut[move_node] = -1;
				sizeCut--;
			}
			//consider its neighborhood
			while (nb < adjacLength[move_node])
			{
				k = adjacTable[move_node][nb];
				if ((set[k]==1) && (degvtoS[k]==0))
				{
					setCut[sizeCut] = k;
					posCut[k] = sizeCut;
					sizeCut++;
				}
				if ((set[k]==0) && (degvtoVS[k]==1))
				{
					setCut[posCut[k]] = setCut[sizeCut-1];
					posCut[setCut[sizeCut-1]] = posCut[k];
					posCut[k] = -1;
					sizeCut--;
				}
				nb++;
			}
		}
	}//end if---from V\S(1) to S(0)

	//from S(0) to V\S(1)
	else
	{
		//not in the cut set
		if (degvtoVS[move_node] == 0)
		{
			//add itself
			setCut[sizeCut] = move_node;
			posCut[move_node] = sizeCut;
			sizeCut++;

			while (nb < adjacLength[move_node])
			{
				k = adjacTable[move_node][nb];
				if (degvtoS[k] == 0)
				{
					fprintf(stderr, "Error in cut_set_update()\n");
					exit(-1);
				}
				if (degvtoVS[k] == 0)
				{
					setCut[sizeCut] = k;
					posCut[k] = sizeCut;
					sizeCut++;
				}
				nb++;
			}
		}
		//in the cut set
		else
		{
			//consider itself
			if (degvtoS[move_node] == 0)
			{
				setCut[posCut[move_node]] = setCut[sizeCut-1];
				posCut[setCut[sizeCut-1]] = posCut[move_node];
				posCut[move_node] = -1;
				sizeCut--;
			}
			//consider its neighborhood
			while (nb < adjacLength[move_node])
			{
				k = adjacTable[move_node][nb];
				if ((set[k]==0) && (degvtoVS[k]==0))
				{
					setCut[sizeCut] = k;
					posCut[k] = sizeCut;
					sizeCut++;
				}
				if ((set[k]==1) && (degvtoS[k]==1))
				{
					setCut[posCut[k]] = setCut[sizeCut-1];
					posCut[setCut[sizeCut-1]] = posCut[k];
					posCut[k] = -1;
					sizeCut--;
				}
				nb++;
			}
		}
	}//end else---from S(0) to V\S(1)
}


/* tabu search */
void cut_neighboring_tabu_search()
{
	int iter;
	long long int add_node;
	long long int eval_degreeS, eval_degreeVS;
	long long int eval_cG, eval_sol;
	long long int best_sol, tabu_best_sol;
	long long int adding_list[10000], tabu_adding_list[10000];
	int nb_best, tabu_nb_best, index;
	int *set;
	long long int sol;

	/* tabu tenure management */
	int p = 0, t = 0;
	const int periods = 15;                      //the number of periods
	const int p_interval = D/periods;           //the intervals: 6000/15=400
	int B[periods] = { 10, 20, 10, 40, 10, 20, 10, 80, 10, 20, 10, 40, 10, 20, 10 };
	int A[periods];                              //alfa
	for (int i=0; i<periods; i++)
		A[i] = alfa;

	//input
	set = new int[nbNodes];
	for (long long int i=0; i<nbNodes; i++)
		set[i] = setFlag[i];
	sol = localBestSol;

	//clear
	iter = 0;
	for (long long int i=0; i<nbNodes; i++)
		tabuTenure[i] = 0;

	//begin
	while (iter < D)
	{
		best_sol = MAX_INT;
		tabu_best_sol = MAX_INT;
		nb_best = 0;
		tabu_nb_best = 0;

		//find(check) the best ADD move in cut edges set
		for (long long int k=0; k<sizeCut; k++)
		{
			//choose the evaluate node
			add_node = setCut[k];

			//calculate the metrics
			if ((set[add_node]==1) && (sizeVS>1))
			{
				eval_degreeS = volumeS + adjacLength[add_node];
				eval_degreeVS = volumeVS - adjacLength[add_node];
				eval_cG = cG - degvtoS[add_node] + degvtoVS[add_node];
			}
			else if ((set[add_node]==0) && (sizeS>1))
			{
				eval_degreeS = volumeS - adjacLength[add_node];
				eval_degreeVS = volumeVS + adjacLength[add_node];
				eval_cG = cG + degvtoS[add_node] - degvtoVS[add_node];
			}
			else
				continue;

			eval_sol = (double)eval_cG/(double)MIN(eval_degreeS, eval_degreeVS)*MAX_INT;

			//if it is tabu
			if (tabuTenure[add_node] > iter)
			{
				if (eval_sol < tabu_best_sol)
				{
					tabu_best_sol = eval_sol;
					tabu_adding_list[0] = add_node;
					tabu_nb_best = 1;
				}
				else if ((eval_sol==tabu_best_sol) && (tabu_nb_best<10000))
				{
					tabu_adding_list[tabu_nb_best] = add_node;
					tabu_nb_best++;
				}
			}
			//if it is not tabu
			else if (tabuTenure[add_node] <= iter)
			{
				if (eval_sol < best_sol)
				{
					best_sol = eval_sol;
					adding_list[0] = add_node;
					nb_best = 1;
				}
				else if ((eval_sol==best_sol) && (nb_best<10000))
				{
					adding_list[nb_best] = add_node;
					nb_best++;
				}
			}
		}

		/*if all nodes in setCut[] are taboo, (special case and rarely happens)
		 * reset tabu tenure to 0, and restart the check
		 */
		if((tabu_nb_best==0) && (nb_best==0))
		{
			for (long long int i=0; i<nbNodes; i++)
				tabuTenure[i] = 0;
			continue;
		}

		//accept tabu solution, i.e. 1)aspiration criterion; 2)no improved sol in non-tabu nodes
		if (((tabu_nb_best>0) && (tabu_best_sol<sol) && (tabu_best_sol<globalBestSol)) || (nb_best==0))
		{
			//randomly select a move
			index = rand()%tabu_nb_best;
			add_node = tabu_adding_list[index];

			sol = tabu_best_sol;

			//make a move
			if (set[add_node] == 1)
			{
				volumeS += adjacLength[add_node];
				volumeVS -= adjacLength[add_node];
				cG = cG - degvtoS[add_node] + degvtoVS[add_node];
				sizeS++;
				sizeVS--;

				cut_set_update(set, add_node);
				set_node_weight(add_node);

				for (long long int i=0; i<adjacLength[add_node]; i++)
				{
					degvtoS[adjacTable[add_node][i]] += 1;
					degvtoVS[adjacTable[add_node][i]] -= 1;
				}

				set[add_node] = 0;
			}
			else
			{
				volumeS -= adjacLength[add_node];
				volumeVS += adjacLength[add_node];
				cG = cG + degvtoS[add_node] - degvtoVS[add_node];
				sizeS--;
				sizeVS++;

				cut_set_update(set, add_node);
				set_node_weight(add_node);

				for (long long int i=0; i<adjacLength[add_node]; i++)
				{
					degvtoS[adjacTable[add_node][i]] -= 1;
					degvtoVS[adjacTable[add_node][i]] += 1;
				}

				set[add_node] = 1;
			}

			//update tabu tenure
			tabuTenure[add_node] = A[p]*B[p];
			tabuTenure[add_node] += iter;
			t++;
			if (t > p_interval)
			{
				p = (p+1) % periods;
				t = 0;
			}
		}

		//accept non-tabu solution (normal condition)
		else
		{
			//randomly select a move
			index = rand()%nb_best;
			add_node = adding_list[index];

			sol = best_sol;

			//make a move
			if (set[add_node] == 1)
			{
				volumeS += adjacLength[add_node];
				volumeVS -= adjacLength[add_node];
				cG = cG - degvtoS[add_node] + degvtoVS[add_node];
				sizeS++;
				sizeVS--;

				cut_set_update(set, add_node);
				set_node_weight(add_node);

				for (long long int i=0; i<adjacLength[add_node]; i++)
				{
					degvtoS[adjacTable[add_node][i]] += 1;
					degvtoVS[adjacTable[add_node][i]] -= 1;
				}

				set[add_node] = 0;
			}
			else
			{
				volumeS -= adjacLength[add_node];
				volumeVS += adjacLength[add_node];
				cG = cG + degvtoS[add_node] - degvtoVS[add_node];
				sizeS--;
				sizeVS++;

				cut_set_update(set, add_node);
				set_node_weight(add_node);

				for (long long int i=0; i<adjacLength[add_node]; i++)
				{
					degvtoS[adjacTable[add_node][i]] -= 1;
					degvtoVS[adjacTable[add_node][i]] += 1;
				}

				set[add_node] = 1;
			}

			//update tabu tenure
			tabuTenure[add_node] = A[p]*B[p];
			tabuTenure[add_node] += iter;
			t++;
			if (t > p_interval)
			{
				p = (p+1) % periods;
				t = 0;
			}
		}

		if (sol < localBestSol)
		{
			localBestSol = sol;
			for (long long int i=0; i<nbNodes; i++)
				setFlag[i] = set[i];
			iter=0;
		}
		else
			iter++;

		if (diff_time(startTime) > runTime)
			break;
	}

	//initialize
	initialize(setFlag, localBestSol);

	//free memory
	delete[] set;
	set = NULL;
}


/* quick sort for weighting nodes */
void quick_sort(int *nod, int *wei, int l, int r)
{
	 if (l < r)
	 {
		int i = l, j = r, x = wei[l], y = nod[l];
		while (i < j)
		{
			//find first < x, from right to left
			while ((i<j) && (wei[j]<=x))
				j--;
			if (i < j)
			{
				wei[i] = wei[j];
				nod[i++] = nod[j];
			}

			//find first >= x, from left to right
			while ((i<j) && (wei[i]>x))
				i++;
			if (i < j)
			{
				wei[j] = wei[i];
				nod[j--] = nod[i];
			}
		}
		wei[i] = x;
		nod[i] = y;
		quick_sort(nod, wei, l, i-1);
		quick_sort(nod, wei, i+1, r);
	}
}


/* type1 (weak): based on ADD operator + vertex weighting scheme */
void local_optimal_exploring_ADD(int &L)
{
	int n, t;
	int *nod, *wei;
	int count;
	long long int add_node;
	int *set;
	long long int sol;

	//randomly select n nodes for weight sorting
	n = L;
	nod = new int[n];
	wei = new int[n];

	for (int i=0; i<n; i++)
	{
		t = (long long int)rand()%nbNodes;
		nod[i] = t;
		wei[i] = weiNodes[t];
	}
	/* quick sort */
	quick_sort(nod, wei, 0, n-1);

	//input
	set = new int[nbNodes];
	for (long long int i=0; i<nbNodes; i++)
		set[i] = setFlag[i];
	sol = localBestSol;

	//begin
	count = 0;
	while (count < L)
	{
		add_node = nod[count];

		//make a move
		if ((set[add_node]==1) && (sizeVS>1))
		{
			volumeS += adjacLength[add_node];
			volumeVS -= adjacLength[add_node];
			cG = cG - degvtoS[add_node] + degvtoVS[add_node];
			sizeS++;
			sizeVS--;

			cut_set_update(set, add_node);
			set_node_weight(add_node);

			for (long long int i=0; i<adjacLength[add_node]; i++)
			{
				degvtoS[adjacTable[add_node][i]] += 1;
				degvtoVS[adjacTable[add_node][i]] -= 1;
			}

			set[add_node] = 0;
		}
		else if ((set[add_node]==0) && (sizeS>1))
		{
			volumeS -= adjacLength[add_node];
			volumeVS += adjacLength[add_node];
			cG = cG + degvtoS[add_node] - degvtoVS[add_node];
			sizeS--;
			sizeVS++;

			cut_set_update(set, add_node);
			set_node_weight(add_node);

			for (long long int i=0; i<adjacLength[add_node]; i++)
			{
				degvtoS[adjacTable[add_node][i]] -= 1;
				degvtoVS[adjacTable[add_node][i]] += 1;
			}

			set[add_node] = 1;
		}
		else
			continue;

		sol = (double)cG/(double)MIN(volumeS, volumeVS)*MAX_INT;

		//record best sol
		if (sol < localBestSol)
		{
			localBestSol = sol;
			for (long long int i=0; i<nbNodes; i++)
				setFlag[i] = set[i];
			count=0;
		}
		else
			count++;

		if (diff_time(startTime) > runTime)
			break;
	}

	//record global best sol
	if (localBestSol < globalBestSol)
	{
		globalBestSol = localBestSol;
		globalTime = diff_time(startTime);
		for (long long int i=0; i<nbNodes; i++)
			globalsetFlag[i] = setFlag[i];
		w=0;
	}
	else
		w++;

	//initialize
	initialize(setFlag, localBestSol);

	//free memory
	delete[] set;
	delete[] nod;
	delete[] wei;
	set = NULL;
	nod = NULL;
	wei = NULL;
}


/* type2 (weak): based on SWAP operator + quality selection method */
void local_optimal_exploring_SWAP(int &L)
{
	int count, num;
	long long int v1, v2;
	long long int eval_degreeS, eval_degreeVS, eval_cG;
	long long int eval_sol, best_sol;
	long long int moving_listS[10000], moving_listVS[10000];
	int nb_best_sol, index;
	int *set;
	long long int sol;

	//input
	set = new int[nbNodes];
	for (long long int i=0; i<nbNodes; i++)
		set[i] = setFlag[i];
	sol = localBestSol;

	//begin
	count = 0;
	while (count < L)
	{
		best_sol= MAX_INT;
		nb_best_sol = 0;

		//find(check) the best SWAP move
		for (v1=0; v1<nbNodes; v1++)
		{
			num = 0;
			while (num < adjacLength[v1])
			{
				v2 = adjacTable[v1][num];
				while ((set[v1]==set[v2]) || (v2<v1))
				{
					num++;
					if(num < adjacLength[v1])
						v2 = adjacTable[v1][num];
					else
						break;
				}

				if ((set[v1]==0) && (set[v2]==1))
				{
					eval_degreeS = volumeS+adjacLength[v2]-adjacLength[v1];
					eval_degreeVS = volumeVS+adjacLength[v1]-adjacLength[v2];
					eval_cG = cG+degvtoS[v1]-degvtoVS[v1]+degvtoVS[v2]-degvtoS[v2]+2;
					eval_sol = (double)eval_cG/(double)MIN(eval_degreeS, eval_degreeVS)*MAX_INT;

					if (eval_sol < best_sol)
					{
						best_sol = eval_sol;
						moving_listS[0] = v1;
						moving_listVS[0] = v2;
						nb_best_sol = 1;
					}
					else if ((eval_sol==best_sol) && (nb_best_sol<10000))
					{
						moving_listS[nb_best_sol] = v1;
						moving_listVS[nb_best_sol] = v2;
						nb_best_sol++;
					}
				}
				else if ((set[v1]==1) && (set[v2]==0))
				{
					eval_degreeS = volumeS+adjacLength[v1]-adjacLength[v2];
					eval_degreeVS = volumeVS+adjacLength[v2]-adjacLength[v1];
					eval_cG = cG+degvtoS[v2]-degvtoVS[v2]+degvtoVS[v1]-degvtoS[v1]+2;
					eval_sol = (double)eval_cG/(double)MIN(eval_degreeS, eval_degreeVS)*MAX_INT;

					if (eval_sol < best_sol)
					{
						best_sol = eval_sol;
						moving_listS[0] = v2;
						moving_listVS[0] = v1;
						nb_best_sol = 1;
					}
					else if ((eval_sol==best_sol) && (nb_best_sol<10000))
					{
						moving_listS[nb_best_sol] = v2;
						moving_listVS[nb_best_sol] = v1;
						nb_best_sol++;
					}
				}

				num++;
			}
		}

		//accept the move without considering the history best sol (***important***)
		if (nb_best_sol > 0)
		{
			//randomly select a move
			index = rand()%nb_best_sol;
			v1 = moving_listS[index];
			v2 = moving_listVS[index];

			//make a move
			volumeS += adjacLength[v2] - adjacLength[v1];
			volumeVS += adjacLength[v1] - adjacLength[v2];
			cG = cG+degvtoS[v1]-degvtoVS[v1]+degvtoVS[v2]-degvtoS[v2]+2;


			//for v1, S to V/S
			cut_set_update(set, v1);
			set_node_weight(v1);
			for (long long int i=0; i<adjacLength[v1]; i++)
			{
				degvtoS[adjacTable[v1][i]] -= 1;
				degvtoVS[adjacTable[v1][i]] += 1;
			}
			set[v1] = 1;

			//for v2, V/S to S
			cut_set_update(set, v2);
			set_node_weight(v2);
			for (long long int i=0; i<adjacLength[v2]; i++)
			{
				degvtoS[adjacTable[v2][i]] += 1;
				degvtoVS[adjacTable[v2][i]] -= 1;
			}
			set[v2] = 0;

			sol = (double)cG/(double)MIN(volumeS, volumeVS)*MAX_INT;
		}

		//record best sol
		if (sol < localBestSol)
		{
			localBestSol = sol;
			for (long long int i=0; i<nbNodes; i++)
				setFlag[i] = set[i];
			count=0;
		}
		else
			count++;

		if (diff_time(startTime) > runTime)
			break;
	}

	//record global best sol
	if (localBestSol < globalBestSol)
	{
		globalBestSol = localBestSol;
		globalTime = diff_time(startTime);
		for (long long int i=0; i<nbNodes; i++)
			globalsetFlag[i] = setFlag[i];
		w=0;
	} else
		w++;

	//initialize
	initialize(setFlag, localBestSol);

	//free memory
	delete[] set;
	set = NULL;
}


/* type3 (strong): based on RANDOM operator */
void local_optimal_escaping(int &L)
{
	int count;
	long long int add_node;

	count = 0;
	while (count < L)
	{
		//randomly select a move
		add_node = (long long int)rand()%nbNodes;

		//make a move
		if ((setFlag[add_node]==1) && (sizeVS>1))
		{
			volumeS += adjacLength[add_node];
			volumeVS -= adjacLength[add_node];
			cG = cG - degvtoS[add_node] + degvtoVS[add_node];
			sizeS++;
			sizeVS--;

			cut_set_update(setFlag, add_node);
			set_node_weight(add_node);

			for (long long int i=0; i<adjacLength[add_node]; i++)
			{
				degvtoS[adjacTable[add_node][i]] += 1;
				degvtoVS[adjacTable[add_node][i]] -= 1;
			}

			setFlag[add_node] = 0;
		}
		else if ((setFlag[add_node]==0) && (sizeS>1))
		{
			volumeS -= adjacLength[add_node];
			volumeVS += adjacLength[add_node];
			cG = cG + degvtoS[add_node] - degvtoVS[add_node];
			sizeS--;
			sizeVS++;

			cut_set_update(setFlag, add_node);
			set_node_weight(add_node);

			for (long long int i=0; i<adjacLength[add_node]; i++)
			{
				degvtoS[adjacTable[add_node][i]] -= 1;
				degvtoVS[adjacTable[add_node][i]] += 1;
			}

			setFlag[add_node] = 1;
		}
		else
			continue;

		localBestSol = (double)cG/(double)MIN(volumeS, volumeVS)*MAX_INT;

		count++;

		if (diff_time(startTime) > runTime)
			break;
	}
}


void SaBTS()
{
//	int count = 0;
	long long int sol;

	//initial and clear
	w = 0;
	L = L0*nbNodes;
	localBestSol = MAX_INT;
	globalBestSol = MAX_INT;
	for (long long int i=0; i<nbNodes; i++)
		weiNodes[i] = 0;

	//record start time
	clock_gettime(CLOCK_REALTIME, &startTime);
	initialize_metis_input();
	cout << "finish initializing" << endl;
	printf("localBestSol = %.10lf\n\n", (double) localBestSol / (double) MAX_INT);

/************************************************************/
	//begin
	while (diff_time(startTime) < runTime)
	{
		/* PART2 improved by tabu search */
		cut_neighboring_tabu_search();

		if (localBestSol < globalBestSol)
		{
			globalBestSol = localBestSol;
			globalTime = diff_time(startTime);
			for (long long int i=0; i<nbNodes; i++)
				globalsetFlag[i] = setFlag[i];
			w=0;
		}
		else
			w++;

/************************************************************/
		/* PART3 determine jump magnitude & perturb types */
		if (w > T)
			w=0;

		if (localBestSol == sol)
			L++;
		else
			L=L0*nbNodes;

		sol = localBestSol;

/************************************************************/
		/* PART4 adaptive diversification strategy */
		if (w == 0)
			local_optimal_escaping(L);
		else
		{
			//calculate p
			double t, p, k, l;
			t = exp(-(double)w/(double)T);
			if ((t-P0) > EPSILON)
				p = t;
			else
				p = P0;

			k = rand_double();

			if ((k-p) < EPSILON)
			{
				l = rand_double();

				if ((l-0.5) < EPSILON)
					local_optimal_exploring_ADD(L);
				else
					local_optimal_exploring_SWAP(L);
			}
			else
				local_optimal_escaping(L);
		}

//		printf("%d - globalBestSol=%.10f, globalTime=%lf\n", count, (double)globalBestSol/(double)MAX_INT, globalTime);
//		count++;

/************************************************************/
        //reset
        for (long long int i=0; i<nbNodes; i++)
        	weiNodes[i] = 0;
	}
}


void free_memory()
{
	delete[] adjacLength;
	delete[] metisFlag;
	delete[] setFlag;
	delete[] globalsetFlag;
	delete[] degvtoS;
	delete[] degvtoVS;
	delete[] setVS;
	delete[] setCut;
	delete[] posCut;
	delete[] weiNodes;
	delete[] tabuTenure;
	adjacLength = NULL;
	metisFlag = NULL;
	setFlag = NULL;
	globalsetFlag = NULL;
	degvtoS = NULL;
	degvtoVS = NULL;
	setVS = NULL;
	setCut = NULL;
	posCut = NULL;
	weiNodes = NULL;
	tabuTenure = NULL;

	for (long long int i = 0; i < nbNodes; i++)
	{
		delete[] adjacTable[i];
		adjacTable[i] = NULL;
	}
	delete[] adjacTable;
	adjacTable = NULL;
}


//for formal run
int main(int argc, char* argv[])
{
	char *instanceName;                   //instance file name
	char *metisName;                      //metis file name
	char *dataSet;                  	  //data set name
	char instancePath[10000];        	  //instance file path
	char metisPath[10000];                // metis file path
	char statisticPath[10000];      	  //statistical file path

	FILE *statistic;
	FILE *out;

	if (argc == 11)
	{
		instanceName = argv[1];
		metisName = argv[2];
		dataSet = argv[3];
		runTime = atoi(argv[4]);
		nbRepeat = atoi(argv[5]);
		alfa = atoi(argv[6]);
		D = atoi(argv[7]);
		T = atoi(argv[8]);
		L0 = atof(argv[9]);
		P0 = atof(argv[10]);
	}
	else
	{
		printf("\n### Input the following parameters ###\n");
		printf("<instanceName> <metisName> <dataSet> <runTime> <nbRepeat> <alfa> <D> <T> <L0> <P0>\n");
		exit(-1);
	}

	//there are 2 types of file path: instances, statistics
	strcpy(instancePath, "./");
	strcat(instancePath, "instances/");
	strcat(instancePath, dataSet);
	strcat(instancePath,"/");
	strcat(instancePath, instanceName);

	strcpy(metisPath, "./");
	strcat(metisPath, "metis/");
	strcat(metisPath, dataSet);
	strcat(metisPath, "/");
	strcat(metisPath, metisName);

	strcpy(statisticPath, "./");
	strcat(statisticPath, "statistics/");
	strcat(statisticPath, dataSet);
	strcat(statisticPath,"/");
	strcat(statisticPath, instanceName);
	strcat(statisticPath, ".statistic");

	if((statistic=fopen(statisticPath, "w")) == NULL)
	{
		printf("Open failed for output %s\n", statisticPath);
		exit(1);
	}

	//read the instance and metis file
	read_graph(instancePath);
	read_metis(metisPath);

	//repeat multiple runs
	long long int rec_sol[nbRepeat];
	double rec_time[nbRepeat];

	fprintf(statistic, "maximum run time=%d(s), number of repeat=%d\n", runTime, nbRepeat);
	fprintf(statistic, "---------------------------------------------------------\n");

	for (int i=0; i<nbRepeat; i++)
	{
		//set random seed
		seed = (unsigned)time(NULL);
		srand(seed);

		//run the SaBTS algorithm
		SaBTS();

		//verify final result
//		check_result();

		//statistical results
		rec_sol[i] = globalBestSol;
		rec_time[i] = globalTime;
		fprintf(statistic, "%.10lf, %.6lf, %6d\n", (double)globalBestSol/(double)MAX_INT, globalTime, i);
	}

	//compute the statistical results
	long long int best_sol=MAX_INT;
	double avg_time=0.0;
	double avg_sol=0.0;
	double sd = 0.0;
	double improved = 0.0;
	int success_count=0;
	long long int worst_sol=0;

	//average data
	for (int i=0; i<nbRepeat; i++)
	{
		avg_time += rec_time[i];
		avg_sol += rec_sol[i];
	}
	avg_time /= nbRepeat;
	avg_sol /= nbRepeat;

	//best and worst data, success rate, sd
	for (int i=0; i<nbRepeat; i++)
	{
		if (rec_sol[i] < best_sol)
			best_sol = rec_sol[i];

		if (rec_sol[i] > worst_sol)
			worst_sol = rec_sol[i];
	}
	for (int i=0; i<nbRepeat; i++)
	{
		if (rec_sol[i] == best_sol)
			success_count++;
	}
	for (int i=0; i<nbRepeat; i++)
		sd += pow((rec_sol[i]-avg_sol), 2);
	sd /= nbRepeat;
	sd = sqrt(sd);
	improved = (double) (metisBestSol - best_sol) / (double) metisBestSol;

	fprintf(statistic,"---------------------------------------------------------\n");
	fprintf(statistic, "%s, %lld, %lld, %lf, %d(s)*%d(times), %.10lf, %.10lf, %.10lf, %.10lf, %lf(s), %d/%d, %.6lf, %.6lf\n",
			instanceName, nbNodes, nbEdges, graphDensity, runTime, nbRepeat,
			(double)metisBestSol/(double)MAX_INT,
			(double)best_sol/(double)MAX_INT,
			(double)worst_sol/(double)MAX_INT,
			(double)avg_sol/(double)MAX_INT, avg_time, success_count, nbRepeat,
			(double)sd/(double)MAX_INT, improved);
	printf("Instance: %s, nodes: %lld, edges: %lld, graph density: %lf\n",
			instanceName, nbNodes, nbEdges, graphDensity);
	printf("The metis sol: %.10f\n", (double) metisBestSol / (double) MAX_INT);
	printf("The best sol: %.10f\n", (double)best_sol/(double)MAX_INT);
	printf("The average sol: %.10f\n", (double)avg_sol/(double)MAX_INT);
	printf("The worst sol: %.10f\n", (double)worst_sol/(double)MAX_INT);
	printf("The hit time: %lf\n", avg_time);
	printf("The success rate: %d/%d\n", success_count, nbRepeat);
	printf("The standard deviation: %.6lf\n", (double)sd/(double)MAX_INT);
	printf("The improved: %.6lf\n", improved);
	fclose(statistic);

	//output for excel
	char finalResultPath[10000];
	strcpy(finalResultPath, "./");
	strcat(finalResultPath, "finalResult.xlsx");
	if ((out=fopen(finalResultPath, "a+")) == NULL)
	{
		printf("Open failed for output %s\n", finalResultPath);
		exit(1);
	}
	fprintf(out, "%s, %lld, %lld, %lf, %d(s)*%d(times), %.10lf, %.10lf, %.10lf, %.10lf, %lf(s), %d//%d, %.6lf, %.6lf\n",
			instanceName, nbNodes, nbEdges, graphDensity, runTime, nbRepeat,
			(double)metisBestSol/(double)MAX_INT,
			(double)best_sol/(double)MAX_INT,
			(double)worst_sol/(double)MAX_INT,
			(double)avg_sol/(double)MAX_INT, avg_time, success_count, nbRepeat,
			(double)sd/(double)MAX_INT, improved);
	fclose(out);

	//free memory
	free_memory();

    return 0;
}


