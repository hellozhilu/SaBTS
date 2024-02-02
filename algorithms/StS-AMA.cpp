//============================================================================
// Name        : StS-AMA.cpp
// Author      : Zhi LU (zhilusix@gmail.com)
// Version     : Feb. 2018
// Copyright   : LERIA, Faculté d'Sciences, Université d'Angers, France
// Description : A Memetic Algorithm for the Minimum Conductance Graph Partitioning Problem.
//               *Designed by David Chalupa (see paper in https://arxiv.org/abs/1704.02854),
//               *Reimplemented by Zhi LU.
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
#define MIN(X,Y)  ((X) <= (Y) ? (X) : (Y))
#define ABS(X)  (((X) < 0) ? -(X) : (X))

using namespace std;


/*graph*/
long long int nbNodes;
long long int nbEdges;
long long int **adjacTable;
long long int *adjacLength;
double graphDensity;

/*global solution*/
long long int volumeS, volumeVS, cG;
long long int sizeS, sizeVS;
long long int localBestSol, globalBestSol;
long long int *degvtoS, *degvtoVS;
long long int *setVS;
int *setFlag, *globalsetFlag;   //labeled vertices belong to which sets, 0-S, 1-V\S
//population
int **popPool;
long long int *popSol;
//crossover
long long int solMom, solDad;
int *mom, *dad;

/*parameters*/
double ps;
int p = 100;
int tournamentSize = 2;
int restartLength = 1;
int maxRLSDepth = 1000000;

/*time*/
int runTime;
int nbRepeat;
struct timespec startTime = {0, 0};
struct timespec endTime = {0, 0};
double globalTime;
long long int seed;


/*function prototype*/
void read_graph(char*);
void allocate_memory();
void check_result();
double diff_time(timespec);
double rand_double();
void restart_set(int*);
void initial_sol(int*, long long int&);
void generate_pop();
void choose_parent();
void crossover();
void move_operator_ADD(int*, long long int&);
void move_operator_SWAP(int*, long long int&);
void LS(int*, long long int&);
void RLS(int*, long long&);
void update_pop();
void StS_AMA();
void free_memory();


void read_graph(char *filename)
{
	string line;
	ifstream ff(filename);
	if (ff.fail())
	{
		fprintf(stderr, "Error in reading file1\n");
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
	long long int *tempLst = new long long int[nbNodes];
	long long int len = 0, lno = 0;
	int vtx;

	while (getline(ff, line))
	{
		if (ff.fail())
		{
			fprintf(stderr, "Error in reading file2\n");
			exit(-1);
		}

		len = 0;
		stringstream ss(line);
		while (ss >> vtx)
			tempLst[len++] = vtx-1;

		adjacLength[lno] = len;
		if (len == 0)
			adjacTable[lno] = NULL;
		else
		{
			adjacTable[lno] = new long long int[len];
			memcpy(adjacTable[lno], tempLst, sizeof(long long int)*len);
		}
		lno++;
	}

	//allocate memory
	allocate_memory();

	//free
	delete[] tempLst;
	tempLst = NULL;

	cout << "finish reading" << endl;
	cout << "nbNodes=" << nbNodes << ", nbEdges=" << nbEdges << ", graphDensity=" << graphDensity << endl << endl;

	ff.close();
}


void allocate_memory()
{
	//initialize other global variables
	setFlag = new int[nbNodes];
	globalsetFlag = new int[nbNodes];
	degvtoS = new long long int[nbNodes];
	degvtoVS = new long long int[nbNodes];
	setVS = new long long int[nbNodes];
	mom = new int[nbNodes];
	dad = new int[nbNodes];
	popSol = new long long int[p];
	popPool = new int*[p];
	for (int i = 0; i < p; i++)
		popPool[i] = new int[nbNodes];
}


void check_result()
{
	long long int min_sol;

	//check
	volumeS=0, volumeVS=0, sizeS=0, sizeVS=0, cG=0;
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

//	spentTime = (double)(temp.tv_sec*1000) + (double)(temp.tv_nsec/1000000);                     //for ms
	spentTime = (double)temp.tv_sec + (double)(temp.tv_nsec/1000)/(double)CLOCKS_PER_SEC;        //for s

	return spentTime;
}


double rand_double()
{
    double r;

    r = (double)rand()/(double)RAND_MAX;
//	r = (rand()%10000)*0.0001;

	return r;
}


void restart_set(int *set)
{
	for (long long int i=0; i<nbNodes; i++)
		set[i] = rand()%2;
}


void initial_sol(int *set, long long int &sol)
{
	//volumeS, volumeVS, cG
	volumeS = 0;
	volumeVS = 0;
	sizeS = 0;
	sizeVS = 0;
	cG = 0;
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
			sizeVS++;
		}
	}

	//ensure the solution is legal
	if ((volumeS==0) || (volumeVS==0))
	{
		printf("volumeS = %lld, volumeVS = %lld\n", volumeS, volumeVS);
		fprintf(stderr, "Error initial solution in initial_sol()\n");
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

	//localBestSol
	sol = (double)cG/(double)MIN(volumeS, volumeVS)*MAX_INT;
}


void generate_pop()
{
	int improve_count;
	double k;
	long long int check_sizeS, check_sizeVS;
	long long int best_sol;

	best_sol = MAX_INT;
	for (int num_pop = 0; num_pop < p; num_pop++)
	{
		ps = (double)1/(double)2;
		improve_count = 0;

		do
		{
			check_sizeS = 0;
			check_sizeVS = 0;
			for (long long int i=0; i<nbNodes; i++)
			{
				k = rand_double();

				if ((k-ps) < EPSILON)
				{
					popPool[num_pop][i] = 1;
					check_sizeVS++;
				}
				else
				{
					popPool[num_pop][i] = 0;
					check_sizeS++;
				}
			}

			//ensure the solution is legal
			if ((check_sizeS==0) || (check_sizeVS==0))
			{
				ps /= 2;
				improve_count++;
				continue;
			}

			initial_sol(popPool[num_pop], popSol[num_pop]);
			/*local search*/
			LS(popPool[num_pop], popSol[num_pop]);

//			printf("num_pop=%d, improve_count=%d\n", num_pop, improve_count);
//			printf("globalBestSol in generate_pop()=%.10f, globalTime=%lf\n\n",
//					(double)globalBestSol/(double)MAX_INT, globalTime);

			ps /= 2;
			improve_count++;

		}while ((popSol[num_pop]>best_sol) && (improve_count<5)
				&& (diff_time(startTime)<runTime));

		//record the previous candidates sol
		if (popSol[num_pop] < best_sol)
			best_sol = popSol[num_pop];
	}
}


void choose_parent()
{
	int num;
	int m, n, t;

	num = 0;
	t = 1000;   //t must > 50
	while (num < tournamentSize)
	{
		m = 0;
		n = 0;
		while ((m==n) || (m==t) || (n==t))
		{
			m = rand()%p;
			n = rand()%p;
		}

		if (num == 0)
		{
			if (popSol[m] > popSol[n])
			{
				for (long long int i=0; i<nbNodes; i++)
					mom[i] = popPool[n][i];
				t = n;
			}
			else
			{
				for (long long int i=0; i<nbNodes; i++)
					mom[i] = popPool[m][i];
				t = m;
			}
		}
		else
		{
			if (popSol[m] > popSol[n])
			{
				for (int i=0; i<nbNodes; i++)
					dad[i] = popPool[n][i];
			}
			else
			{
				for (int i=0; i<nbNodes; i++)
					dad[i] = popPool[m][i];
			}
		}
		num++;
	}
}


//one-point crossover
void crossover()
{
	int pos;
	int *part;

	//except 0, nbNodes
	pos = 0;
	while (pos == 0)
		pos = rand()%nbNodes;

	part = new int[pos];

	for (int i=0; i<pos; i++)
	{
		part[i] = mom[i];
		mom[i] = dad[i];
		dad[i] = part[i];
	}

	//free
	delete[] part;
	part = NULL;
}


void move_operator_ADD(int *set, long long int &sol)
{
	long long int total_degreeS, total_degreeVS;
	long long int eval_cG, eval_sol;
	long long int adding_list[10000];
	long long int best_sol, move_node;
	int nb_best_sol, index;

	while (1)
	{
		best_sol = MAX_INT;
		nb_best_sol = 0;

		//find the best neighboring move
		for (long long int k=0; k<nbNodes; k++)
		{
			//calculate the metrics
			if ((set[k]==1) && (sizeVS>1))
			{
				total_degreeS = volumeS + adjacLength[k];
				total_degreeVS = volumeVS - adjacLength[k];
				eval_cG = cG - degvtoS[k] + degvtoVS[k];
			}
			else if ((set[k]==0) && (sizeS>1))
			{
				total_degreeS = volumeS - adjacLength[k];
				total_degreeVS = volumeVS + adjacLength[k];
				eval_cG = cG + degvtoS[k] - degvtoVS[k];
			}
			else
				continue;

			eval_sol = (double)eval_cG/(double)MIN(total_degreeS, total_degreeVS)*MAX_INT;

			if (eval_sol < best_sol)
			{
				best_sol = eval_sol;
				adding_list[0] = k;
				nb_best_sol = 1;
			}
			else if ((eval_sol==best_sol) && (nb_best_sol<10000))
			{
				adding_list[nb_best_sol] = k;
				nb_best_sol++;
			}
		}

		//update solution
		if ((best_sol<sol) && (nb_best_sol>0))
		{
			//randomly select a move
			index = rand()%nb_best_sol;
			move_node = adding_list[index];

			sol = best_sol;

			//make a move
			if (set[move_node] == 1)
			{
				volumeS += adjacLength[move_node];
				volumeVS -= adjacLength[move_node];
				cG = cG - degvtoS[move_node] + degvtoVS[move_node];
				sizeS++;
				sizeVS--;

				for (long long int i=0; i<adjacLength[move_node]; i++)
				{
					degvtoS[adjacTable[move_node][i]] += 1;
					degvtoVS[adjacTable[move_node][i]] -= 1;
				}

				set[move_node] = 0;
			}
			else
			{
				volumeS -= adjacLength[move_node];
				volumeVS += adjacLength[move_node];
				cG = cG + degvtoS[move_node] - degvtoVS[move_node];
				sizeS--;
				sizeVS++;

				for (long long int i=0; i<adjacLength[move_node]; i++)
				{
					degvtoS[adjacTable[move_node][i]] -= 1;
					degvtoVS[adjacTable[move_node][i]] += 1;
				}

				set[move_node] = 1;
			}
		}
		else
			break;

		if (diff_time(startTime) > runTime)
			break;
	}
}

//new
void move_operator_SWAP(int *set, long long int &sol)
{
	long long int num, v0, v1;
	long long int eval_degreeS, eval_degreeVS, eval_cG;
	long long int eval_sol, best_sol;
	long long int moving_listS[10000], moving_listVS[10000];
	int nb_best_sol, index;

	//begin
	while (1)
	{
		best_sol= MAX_INT;
		nb_best_sol = 0;

		for (v0=0; v0<nbNodes; v0++)
		{
			num = 0;
			while (num < adjacLength[v0])
			{
				v1 = adjacTable[v0][num];
				while ((set[v0]==set[v1]) || (v1<v0))
				{
					num++;
					if(num < adjacLength[v0])
						v1 = adjacTable[v0][num];
					else
						break;
				}

				if ((set[v0]==0) && (set[v1]==1))
				{
					eval_degreeS = volumeS+adjacLength[v1]-adjacLength[v0];
					eval_degreeVS = volumeVS+adjacLength[v0]-adjacLength[v1];
					eval_cG = cG+degvtoS[v0]-degvtoVS[v0]+degvtoVS[v1]-degvtoS[v1]+2;
					eval_sol = (double)eval_cG/(double)MIN(eval_degreeS, eval_degreeVS)*MAX_INT;

					if (eval_sol < best_sol)
					{
						best_sol = eval_sol;
						moving_listS[0] = v0;
						moving_listVS[0] = v1;
						nb_best_sol = 1;
					}
					else if ((eval_sol==best_sol) && (nb_best_sol<10000))
					{
						moving_listS[nb_best_sol] = v0;
						moving_listVS[nb_best_sol] = v1;
						nb_best_sol++;
					}
				}
				else if ((set[v0]==1) && (set[v1]==0))
				{
					eval_degreeS = volumeS+adjacLength[v0]-adjacLength[v1];
					eval_degreeVS = volumeVS+adjacLength[v1]-adjacLength[v0];
					eval_cG = cG+degvtoS[v1]-degvtoVS[v1]+degvtoVS[v0]-degvtoS[v0]+2;
					eval_sol = (double)eval_cG/(double)MIN(eval_degreeS, eval_degreeVS)*MAX_INT;

					if (eval_sol < best_sol)
					{
						best_sol = eval_sol;
						moving_listS[0] = v1;
						moving_listVS[0] = v0;
						nb_best_sol = 1;
					}
					else if ((eval_sol==best_sol) && (nb_best_sol<10000))
					{
						moving_listS[nb_best_sol] = v1;
						moving_listVS[nb_best_sol] = v0;
						nb_best_sol++;
					}
				}

				num++;
			}
		}

		//make moves
		if ((best_sol<sol) && (nb_best_sol>0))
		{
			//randomly select moves
			index = rand()%nb_best_sol;
			v0 = moving_listS[index];
			v1 = moving_listVS[index];

			//make moves
			volumeS += adjacLength[v1] - adjacLength[v0];
			volumeVS += adjacLength[v0] - adjacLength[v1];
			cG = cG+degvtoS[v0]-degvtoVS[v0]+degvtoVS[v1]-degvtoS[v1]+2;

			//for v0, S to V/S
			for (long long int i=0; i<adjacLength[v0]; i++)
			{
				degvtoS[adjacTable[v0][i]] -= 1;
				degvtoVS[adjacTable[v0][i]] += 1;
			}
			set[v0] = 1;
			//for v1, V/S to S
			for (long long int i=0; i<adjacLength[v1]; i++)
			{
				degvtoS[adjacTable[v1][i]] += 1;
				degvtoVS[adjacTable[v1][i]] -= 1;
			}
			set[v1] = 0;

			sol = (double)cG/(double)MIN(volumeS, volumeVS)*MAX_INT;
		}
		else
			break;

		if (diff_time(startTime) > runTime)
			break;
	}
}


void LS(int *setFlag, long long int &localBestSol)
{
	int count;
	int *set;
	long long int sol;

	set = new int[nbNodes];
	for (long long int i=0; i<nbNodes; i++)
		set[i] = setFlag[i];
	sol = localBestSol;

	count = 0;
	while (count < restartLength)
	{
		move_operator_ADD(set, sol);

		if (sol < localBestSol)
		{
			localBestSol = sol;
			for (long long int i=0; i<nbNodes; i++)
				setFlag[i] = set[i];
			count = 0;
		}
		else
			count++;

		//restart
		restart_set(set);
		initial_sol(set, sol);

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
	}

	//free
	delete [] set;
	set = NULL;
}


void RLS(int *setFlag, long long int &localBestSol)
{
	int no_improved;
	int *set;
	long long int sol;

	set = new int[nbNodes];
	for (long long int i=0; i<nbNodes; i++)
		set[i] = setFlag[i];
	sol = localBestSol;

	for (int count = 0; count < restartLength; count++)
	{
		no_improved = 0;
		while (no_improved < maxRLSDepth)
		{

			move_operator_ADD(set, sol);

			move_operator_SWAP(set, sol);

			if (sol < localBestSol)
			{
				localBestSol = sol;
				for (long long int i=0; i<nbNodes; i++)
					setFlag[i]= set[i];
				no_improved = 0;
			}
			else
				no_improved++;

			if (diff_time(startTime) > runTime)
				break;
		}

		//restart
		restart_set(set);
		initial_sol(set, sol);

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
	}

	//free
	delete [] set;
	set = NULL;
}


void update_pop()
{
	int flag0, flag1, index;
	long long int remove_node, worst_sol, best_sol;
	int worst_nb, worst_index[10000];

	//O1, O2 must NOT belong to the population
	for (int i=0; i<p; i++)
	{
		if (popSol[i] != solMom)
			flag0 = 0;
		else
			flag0 = 1;

		if (popSol[i] != solDad)
			flag1 = 0;
		else
			flag1 = 1;
	}

	//replace the worst individual in P with O1
	if (flag0 == 0)
	{
		worst_sol = 0;
		for (int i=0; i<p; i++)
		{
			if (popSol[i] > worst_sol)
			{
				worst_sol = popSol[i];
				worst_index[0] = i;
				worst_nb = 1;
			}
			else if ((popSol[i]==worst_sol) && (worst_nb<10000))
			{
				worst_index[worst_nb] = i;
				worst_nb++;
			}
		}

		if (solMom < worst_sol)
		{
			index = rand()%worst_nb;
			remove_node = worst_index[index];

			for (int i=0; i<nbNodes; i++)
				popPool[remove_node][i] = mom[i];
			popSol[remove_node] = solMom;
		}
	}
	//replace the worst individual in P with O2
	if (flag1 == 0)
	{
		worst_sol = 0;
		for (int i=0; i<p; i++)
		{
			if (popSol[i] > worst_sol)
			{
				worst_sol = popSol[i];
				worst_index[0] = i;
				worst_nb = 1;
			}
			else if ((popSol[i]==worst_sol) && (worst_nb<10000))
			{
				worst_index[worst_nb] = i;
				worst_nb++;
			}
		}

		if (solDad < worst_sol)
		{
			index = rand()%worst_nb;
			remove_node = worst_index[index];

			for (int i=0; i<nbNodes; i++)
				popPool[remove_node][i] = dad[i];
			popSol[remove_node] = solDad;
		}
	}

    //return the best individual in population
	best_sol = MAX_INT;
	for (int i=0; i<p; i++)
	{
		if (popSol[i] < best_sol)
		{
			best_sol = popSol[i];
			index = i;
		}
	}

	for (int i=0; i<nbNodes; i++)
		setFlag[i] = popPool[index][i];
	localBestSol = popSol[index];
}


void StS_AMA()
{
	int gen = 0;
	globalBestSol = MAX_INT;

	//record start time
	clock_gettime(CLOCK_REALTIME, &startTime);

	/*0. generate initial population and improved using LS*/
	generate_pop();
	cout << "finish initializing" << endl;
	printf("globalBestSol=%.10f, globalTime=%lf\n\n", (double)globalBestSol/(double)MAX_INT, globalTime);

	while (diff_time(startTime) < runTime)
	{
		choose_parent();

		//ensure the solution is legal
		int sizeS0 = 0, sizeS1 = 0;
		while((sizeS0==0||sizeS0==nbNodes) || (sizeS1==0||sizeS1==nbNodes))
		{
			crossover();
			sizeS0 = 0, sizeS1 = 0;
			for(long long int i=0; i<nbNodes; i++)
			{
				if(mom[i] == 0)
					sizeS0++;
				if(dad[i] == 0)
					sizeS1++;
			}
		}

		/*1. improve the offspring O1, O2 using ARLS for l = 1000000 iterations*/
		initial_sol(mom, solMom);  //for mom[]
		RLS(mom, solMom);
		initial_sol(dad, solDad);  //for dad[]
		RLS(dad, solDad);

		/*2. improve the offspring O1, O2 using LS until the local optimum is reached*/
		initial_sol(mom, solMom);  //for mom[]
		LS(mom, solMom);
		initial_sol(dad, solDad);  //for dad[]
		LS(dad, solDad);

		/*3. update population*/
		update_pop();

		if (localBestSol < globalBestSol)
		{
			globalBestSol = localBestSol;
			globalTime = diff_time(startTime);;
			for (long long int i=0; i<nbNodes; i++)
				globalsetFlag[i] = setFlag[i];
		}

//		printf("%d - globalBestSol=%.10f, globalTime=%lf\n", gen, (double)globalBestSol/(double)MAX_INT, globalTime);
		gen++;
	}
}


void free_memory()
{
	delete[] adjacLength;
	delete[] degvtoS;
	delete[] degvtoVS;
	delete[] setVS;
	delete[] mom;
	delete[] dad;
	delete[] setFlag;
	delete[] globalsetFlag;
	delete[] popSol;
	adjacLength = NULL;
	degvtoS = NULL;
	degvtoVS = NULL;
	setVS = NULL;
	mom = NULL;
	dad = NULL;
	setFlag = NULL;
	globalsetFlag = NULL;
	popSol = NULL;
	for (int i = 0; i < p; i++)
	{
		delete popPool[i];
		popPool[i] = NULL;
	}
	delete[] popPool;
	popPool = NULL;
	for (long long int i = 0; i < nbNodes; i++)
	{
		delete[] adjacTable[i];
		adjacTable[i] = NULL;
	}
	delete[] adjacTable;
	adjacTable = NULL;
}


int main(int argc, char* argv[])
{
	char *instanceName; 		    //instance file name
	char *dataSet;                  //data set name
	char instancePath[10000];		//instance file path
	char statisticPath[10000];		//statistical file path

	FILE *statistic;
	FILE *out;

	if (argc == 5)
	{
		instanceName = argv[1];
		dataSet = argv[2];
		runTime = atoi(argv[3]);
		nbRepeat = atoi(argv[4]);
	}
	else
	{
		printf("\n### Input the following parameters ###\n");
		printf("<instanceName> <dataSet> <runTime> <nbRepeat>\n");
		exit(-1);
	}

	/* there are 2 types of file path: instances, statistics */
	strcpy(instancePath, "./");
	strcat(instancePath, "instances/");
	strcat(instancePath, dataSet);
	strcat(instancePath,"/");
	strcat(instancePath, instanceName);

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

	//read the instance
	read_graph(instancePath);

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

		//run the StS_AMA algorithm
		StS_AMA();

		//verify final result
//		check_result();

		//statistical results
		rec_sol[i] = globalBestSol;
		rec_time[i] = globalTime;
		fprintf(statistic, "%.10lf, %.6lf, %6d\n", (double)globalBestSol/(double)MAX_INT, globalTime, i);
	}

	/* compute the statistical results */
	long long int best_sol=MAX_INT;
	double avg_time=0.0;
	double avg_sol=0.0;
	double sd = 0.0;
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

	fprintf(statistic,"---------------------------------------------------------\n");
	fprintf(statistic, "%s, %lld, %lld, %lf, %d(s)*%d(times), %.10lf, %.10lf, %.10lf, %lf(s), %d/%d, %lf\n",
			instanceName, nbNodes, nbEdges, graphDensity, runTime, nbRepeat,
			(double)best_sol/(double)MAX_INT,
			(double)worst_sol/(double)MAX_INT,
			(double)avg_sol/(double)MAX_INT, avg_time, success_count, nbRepeat,
			(double)sd/(double)MAX_INT);
	printf("Instance: %s, nodes: %lld, edges: %lld, graph density: %lf\n",
			instanceName, nbNodes, nbEdges, graphDensity);
	printf("The best sol: %.10f\n", (double)best_sol/(double)MAX_INT);
	printf("The average sol: %.10f\n", (double)avg_sol/(double)MAX_INT);
	printf("The worst sol: %.10f\n", (double)worst_sol/(double)MAX_INT);
	printf("The hit time: %lf\n", avg_time);
	printf("The success rate: %d/%d\n", success_count, nbRepeat);
	printf("The standard deviation: %lf\n", (double)sd/(double)MAX_INT);
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
	fprintf(out, "%s, %lld, %lld, %lf, %d(s)*%d(times), %.10lf, %.10lf, %.10lf, %lf(s), %d/%d, %lf\n",
			instanceName, nbNodes, nbEdges, graphDensity, runTime, nbRepeat,
			(double)best_sol/(double)MAX_INT,
			(double)worst_sol/(double)MAX_INT,
			(double)avg_sol/(double)MAX_INT, avg_time, success_count, nbRepeat,
			(double)sd/(double)MAX_INT);
	fclose(out);

	//free memory
	free_memory();

    return 0;
}
