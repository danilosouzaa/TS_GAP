/*
 ============================================================================
 name        : TS-GAP.cu
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : CUDA compute reciprocals
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <curand.h>
#include <curand_kernel.h>
#include <sys/time.h>


#include "gpulib/types.h"
#include "gpulib/gpu.cuh"

#include "Instance.h"
#include "Solution.h"
#include "gSolution.cuh"
#include "guloso.h"

const int nThreads = 100;
const int nBlocks = 3;
const int maxChain = 10;

int main(int argc, char *argv[])
{
	//Variable with GPU's number
	int deviceCount = 0;
	//Commands for verify use correct of GPU
	cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
	if (error_id != cudaSuccess)
	{
		printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
		printf("Result = FAIL\n");
		exit(1);
	}
	if(deviceCount == 0)
	{
		printf("No GPU found :(");
		exit(1);
	}
	else
	{
		printf("Found %d GPUs!\n", deviceCount);
		gpuSetDevice(0);
		printf("GPU 0 initialized!\n");
	}
	//iterator of use in for
	int i, j;
	//Pointer of states for use in curand (GPU)
	curandState_t *states;
	cudaMalloc((void**)&states, (nThreads*nBlocks)*sizeof(curandState_t));

	//Pointer of seed for use with curand (host)
	unsigned int *h_seed = (unsigned int*)malloc(sizeof(unsigned int)*(nThreads*nBlocks));
	srand(time(NULL));
	for(i=0;i<(nThreads*nBlocks);i++){
		h_seed[i] = rand()%100000;
	}


	//Pointers of intance, solution and ejection for use in CPU(Host) and GPU(Device)
	Instance *h_instance, *d_instance;
	Solution *h_solution, *d_solution;
	EjectionChain *h_ejection, *d_ejection;

	//Load Instance
	char nameAux[50] = "../Instances/";
	const char *temp_teste = argv[1];
	strcat(nameAux,argv[1]);
	printf("%s \n",nameAux);
	const char *fileName = nameAux; //argv[1];
//	strcat(fileName,argv[1]);
	h_instance = loadInstance(fileName);

	//showInstance(h_instance);
	//Allocation Solution and Ejection
	h_solution = allocationPointersSolution(h_instance);
	h_ejection = allocationPointerEjectionChain(h_instance);
	//weight greedy
	float w1,w2;
	struct timeval time_rand;

	//Generate Initial Solution from greedy method
	for(i=0;i<nBlocks;i++){
		gettimeofday(&time_rand,NULL);
		srand(time_rand.tv_usec);
		//memset(h_solution,0,size_solution);

		if(temp_teste[0]=='e'){
			do{
				for(j=0;j<h_instance->mAgents;j++){
					h_solution->resUsage[j+i*h_instance->mAgents] = 0;
				}
				w1 = (float)(rand())/(float)(RAND_MAX) + 0.5;
				w2 = 19 + w1;
			}while(guloso(h_instance,h_solution,w1,w2,i)==0);
		}else{
			do{
				for(j=0;j<h_instance->mAgents;j++){
					h_solution->resUsage[j+i*h_instance->mAgents] = 0;
				}
				w1 = (float)(rand())/(float)(RAND_MAX) + 0.5;
				w2 = 1 + w1;
			}while(guloso(h_instance,h_solution,w1,w2,i)==0);
		}
	}
	//Size Struct Solution
	size_t size_solution =  sizeof(Solution)
																			+ sizeof(TcostFinal)*nBlocks
																			+ sizeof(Ts)*(h_instance->nJobs*nBlocks)
																			+ sizeof(TresUsage)*(h_instance->mAgents*nBlocks);

	//Size Struct of Ejection Chain
	size_t size_ejection = sizeof(EjectionChain)
											+ sizeof(Tpos)*(nBlocks*nThreads*maxChain)
											+ sizeof(Top)*(nBlocks*nThreads)
											+ sizeof(TSizeChain)*(nBlocks*nThreads)
											+ sizeof(Tdelta)*(nBlocks*nThreads);

	//Size Struct of Instance
	size_t size_instance = sizeof(Instance)
																									+ sizeof(Tcost)*(h_instance->nJobs*h_instance->mAgents)  //cost
																									+ sizeof(TresourcesAgent)*(h_instance->nJobs*h_instance->mAgents)
																									+ sizeof(Tcapacity)*h_instance->mAgents;

	int *h_short_list = (int*)malloc(sizeof(int)*(nBlocks*h_instance->nJobs));
	int *h_long_list = (int*)malloc(sizeof(int)*(h_instance->nJobs*h_instance->mAgents));
	memset(h_short_list,0,sizeof(int)*(nBlocks*h_instance->nJobs));
	memset(h_long_list,0,sizeof(int)*(h_instance->nJobs*h_instance->mAgents));
	int cost_saida = 1000000;
	for(i=0;i<nBlocks;i++){
		for(j=0;j<h_instance->nJobs;j++){
			h_long_list[j + h_solution->s[j+i*h_instance->nJobs]*h_instance->nJobs]++;
		}
	}

	for(i=0;i<nBlocks;i++){
		printf("Initial cost: %d\n", h_solution->costFinal[i]);
		if(cost_saida>h_solution->costFinal[i]){
			cost_saida=h_solution->costFinal[i];
		}
//		if(i==0){
//			for(j=0;j<h_instance->nJobs;j++){
//				printf("job %d agent %d\n",j, h_solution->s[j+i*h_instance->nJobs]);
//			}			
//		}
	}

	int nJ = h_instance->nJobs;

	int *d_short_list;
	gpuMalloc((void*)&d_short_list,sizeof(int)*(nBlocks*h_instance->nJobs) );
	gpuMemcpy(d_short_list, h_short_list,sizeof(int)*(nBlocks*h_instance->nJobs), cudaMemcpyHostToDevice);


	int blockSize;      // The launch configurator returned block size
	int minGridSize;    // The minimum grid size needed to achieve the maximum occupancy for a full device launch
	int gridSize;
	int N = 1000000;

	cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize,TS_GAP, 0, N);

		printf("block size %d\n",blockSize);
		printf("Min Grid %d\n",minGridSize);
	//	getchar();



	//Reallocation of pointers Instance and Solution for GPU (device)
	d_instance = createGPUInstance(h_instance, h_instance->nJobs, h_instance->mAgents);
	d_solution = createGPUsolution(h_solution,h_instance->nJobs, h_instance->mAgents);
	d_ejection = createGPUejection(h_ejection,h_instance->nJobs, h_instance->mAgents);

	//Pointers seed in device (GPU)
	unsigned int *d_seed;

	//Event and gpu for contability time
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	// Allocation of pointer and copy value in d_seed (Device)
	gpuMalloc((void*)&d_seed, sizeof(unsigned int)*(nThreads*nBlocks));
	gpuMemcpy(d_seed, h_seed, sizeof(unsigned int)*(nThreads*nBlocks), cudaMemcpyHostToDevice);

	cudaEventRecord(start);
	int n_iteration = atoi(argv[2]);
	int ite=1;
	int n_busca = atoi(argv[3]);
	int b_solution = atoi(argv[4]);
	int ite_b = 0;
	int sizeTabu;
	int menor,aux1,t1,m1,m2,aux;
	int *v_menor_pos = (int*)malloc(sizeof(int)*nBlocks);
	
	struct timeval inicio;
	struct timeval fim;
	int tmili = 0;

	nJ = nJ/(maxChain+1);
	gettimeofday(&inicio, NULL);
	while((ite<=n_iteration)&&(tmili<=1200000)){
		sizeTabu = rand()%nJ + 1;
		printf("Size tabu: %d\n", sizeTabu);
		TS_GAP<<<nBlocks,nThreads>>>(d_instance, d_solution,d_ejection, d_short_list, d_seed, states, ite, n_busca);
		cudaDeviceSynchronize();
		gpuMemcpy(h_instance, d_instance, size_instance, cudaMemcpyDeviceToHost);
		gpuMemcpy(h_solution, d_solution, size_solution, cudaMemcpyDeviceToHost);
		gpuMemcpy(h_ejection, d_ejection, size_ejection, cudaMemcpyDeviceToHost);
		gpuMemcpy(h_short_list, d_short_list,sizeof(int)*(nBlocks*h_instance->nJobs), cudaMemcpyDeviceToHost);
		gpuMemcpy(h_seed, d_seed, sizeof(unsigned int)*(nThreads*nBlocks), cudaMemcpyDeviceToHost);

		//reallocation pointers of Instance
		h_instance->cost = (Tcost*)(h_instance+1);
		h_instance->resourcesAgent =(TresourcesAgent*) (h_instance->cost +(h_instance->nJobs*h_instance->mAgents));
		h_instance->capacity =(Tcapacity*) (h_instance->resourcesAgent + (h_instance->nJobs*h_instance->mAgents));

		//reallocation pointers of Solution
		h_solution->costFinal = (TcostFinal*)(h_solution+1);
		h_solution->s = (Ts*)(h_solution->costFinal + nBlocks);
		h_solution->resUsage = (TresUsage*)(h_solution->s + (h_instance->nJobs*nBlocks));


		//reallocation pointers of Ejection
		h_ejection->pos=(Tpos*)(h_ejection + 1);
		h_ejection->op = (Top*)(h_ejection->pos+ (nBlocks*nThreads*maxChain));
		h_ejection->sizeChain = (TSizeChain*)(h_ejection->op + (nBlocks*nThreads));
		h_ejection->delta = (Tdelta*)(h_ejection->sizeChain + (nBlocks*nThreads));


		printf("%d time %d \n",ite,tmili);
		for(i=0;i<nBlocks;i++){
			menor = 100000;
			for(j=0;j<nThreads;j++){
				if(h_ejection->delta[j + i*nThreads]<menor){
					menor = h_ejection->delta[j + i*nThreads];
				}
				//printf("value of delta for thread %d in block %d: :%d \n", j, i, h_ejection->delta[j + i*nThreads]);
			}
			menor = returnIndice(h_solution,h_ejection,i,nBlocks,nThreads,menor,h_long_list,h_instance->nJobs,h_instance->mAgents);
		//	printf("menor delta do bloco %d: %d\n",i,menor);
			if(h_ejection->op[menor + i*nThreads]==1){
				aux1 = h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads];
				//aux2 = ejection->pos[0 + menor*maxChain + i*maxChain*nThreads];
				h_short_list[aux1 + i*h_instance->nJobs] = ite + sizeTabu;

			}else{
				for(j = 0; j<h_ejection->sizeChain[menor + i*nThreads];j++){
					aux1 = h_ejection->pos[j + menor*maxChain + i*maxChain*nThreads];
					h_short_list[aux1 + i*h_instance->nJobs] = ite + sizeTabu;
				}

			}

			h_solution->costFinal[i] += h_ejection->delta[menor+i*nThreads];
			if(h_ejection->op[menor + i*nThreads]==1){
				t1 = h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads];
				m2 = h_ejection->pos[1 + menor*maxChain + i*maxChain*nThreads];
				m1 = ((int)h_solution->s[t1 + i*h_instance->nJobs]);
				h_solution->resUsage[m1 + i*h_instance->mAgents] -= h_instance->resourcesAgent[t1*h_instance->mAgents + m1];
				h_solution->resUsage[m2 + i*h_instance->mAgents] += h_instance->resourcesAgent[t1*h_instance->mAgents + m2];
				h_solution->s[t1 + i*h_instance->nJobs] = ((char)m2);
		//		if(m2>4){
		//			printf("op 1");
		//		}

			}else{
				h_solution->resUsage[((int)h_solution->s[h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs]) + i*h_instance->mAgents] += h_instance->resourcesAgent[h_ejection->pos[(h_ejection->sizeChain[menor + i*nThreads]-1)  + menor*maxChain + i*maxChain*nThreads]*h_instance->mAgents + ((int)h_solution->s[h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs])];
				h_solution->resUsage[((int)h_solution->s[h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs]) + i*h_instance->mAgents] -= h_instance->resourcesAgent[h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads]*h_instance->mAgents + ((int)h_solution->s[h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs])];
				aux = ((int)h_solution->s[h_ejection->pos[0 + menor*maxChain + i*maxChain*nThreads]+ i*h_instance->nJobs]);
				for(j=1; j<h_ejection->sizeChain[menor + i*nThreads]; j++){
					h_solution->resUsage[((int)h_solution->s[h_ejection->pos[j + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs]) + i*h_instance->mAgents] += h_instance->resourcesAgent[h_ejection->pos[(j-1) + menor*maxChain + i*maxChain*nThreads]*h_instance->mAgents + ((int)h_solution->s[h_ejection->pos[j + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs])];
					h_solution->resUsage[((int)h_solution->s[h_ejection->pos[j + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs]) + i*h_instance->mAgents] -= h_instance->resourcesAgent[h_ejection->pos[j + menor*maxChain + i*maxChain*nThreads]*h_instance->mAgents + ((int)h_solution->s[h_ejection->pos[j + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs])];
					h_solution->s[h_ejection->pos[(j-1) + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs] = h_solution->s[h_ejection->pos[j + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs];
				}
				h_solution->s[h_ejection->pos[(h_ejection->sizeChain[menor + i*nThreads]-1) + menor*maxChain + i*maxChain*nThreads] + i*h_instance->nJobs] = ((char)aux);
			}
//			printf("cost: %d\n", h_solution->costFinal[i]);
			      if(cost_saida>h_solution->costFinal[i]){
		                        cost_saida = h_solution->costFinal[i];
                		}

		}

		for(i=0;i<nBlocks;i++){
			for(j=0;j<h_instance->nJobs;j++){
				h_long_list[j + h_solution->s[j+i*h_instance->nJobs]*h_instance->nJobs]++;
			}
		}


	/*	for(i=0;i<nBlocks;i++){
			for(j=0;j<h_instance->nJobs;j++){
				h_long_list[j + h_solution->s[j+i*h_instance->nJobs]*h_instance->nJobs]++;
				if(h_solution->s[j + i*h_instance->nJobs]>4){
					printf("cpu teste: %d\n",h_solution->s[j + i*h_instance->nJobs]);
				}
			}
		}*/


		if(ite!=n_iteration){
			//reallocation pointers of Instance
			h_instance->cost = (Tcost*)(d_instance+1);
			h_instance->resourcesAgent =(TresourcesAgent*) (h_instance->cost +(h_instance->nJobs*h_instance->mAgents));
			h_instance->capacity =(Tcapacity*) (h_instance->resourcesAgent + (h_instance->nJobs*h_instance->mAgents));
			gpuMemcpy(d_instance, h_instance,size_instance, cudaMemcpyHostToDevice);

			//reallocation pointers of Solution
			h_solution->costFinal = (TcostFinal*)(d_solution+1);
			h_solution->s = (Ts*)(h_solution->costFinal + nBlocks);
			h_solution->resUsage = (TresUsage*)(h_solution->s + (h_instance->nJobs*nBlocks));
			gpuMemcpy(d_solution, h_solution, size_solution, cudaMemcpyHostToDevice);

			//reallocation pointers of Ejection
			memset(h_ejection,0,size_ejection);
			h_ejection->pos=(Tpos*)(d_ejection + 1);
			h_ejection->op = (Top*)(h_ejection->pos+ (nBlocks*nThreads*maxChain));
			h_ejection->sizeChain = (TSizeChain*)(h_ejection->op + (nBlocks*nThreads));
			h_ejection->delta = (Tdelta*)(h_ejection->sizeChain + (nBlocks*nThreads));
			gpuMemcpy(d_ejection, h_ejection, size_ejection, cudaMemcpyHostToDevice);

			gettimeofday(&time_rand,NULL);
			srand(time_rand.tv_usec);
			for(i=0;i<(nThreads*nBlocks);i++){
				h_seed[i] = rand()%100000;
			}
			gpuMemcpy(d_seed, h_seed, sizeof(unsigned int)*(nThreads*nBlocks), cudaMemcpyHostToDevice);
			gpuMemcpy(d_short_list, h_short_list,sizeof(int)*(nBlocks*h_instance->nJobs), cudaMemcpyHostToDevice);
			
		}
		if((ite_b==0)&&(cost_saida<= 1.01 * b_solution)){
			ite_b = ite;
		}
		ite++;
		gettimeofday(&fim, NULL);
		tmili = (int) (1000 * (fim.tv_sec - inicio.tv_sec) + (fim.tv_usec - inicio.tv_usec) / 1000);
	}
	printf("cost: %d\n ite: %d\n", cost_saida, ite_b);

	cudaFree(states);
	cudaFree(d_instance);
	cudaFree(d_solution);
	cudaFree(d_ejection);
	cudaFree(d_seed);
	cudaFree(d_short_list);

	free(h_short_list);
	free(h_seed);
	free(h_instance);
	free(h_solution);
	free(h_ejection);
	return 0;
}


