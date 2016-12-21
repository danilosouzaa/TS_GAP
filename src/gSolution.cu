#include "gSolution.cuh"

const int nBlocks = 24;
const int nThreads = 576;
const int maxChain = 8;


__global__ void TS_GAP(Instance *inst, Solution *sol,EjectionChain *ejection, int *tabuListshort, unsigned int *seed, curandState_t* states, int iteration,int sizeTabu, int n_busca)
{
	int i,j,k,flag,a,w;
	//variable for verify if solution is feasible
	int aux = 0;
	int res_aux[80];
	int delta_aux[maxChain];
	int size_aux;
	int d_best_aux;
	int delta_best;
	short int pos_best[maxChain];
	short int size_best;
	short int op_best;

	//initializes of curand
	curand_init(seed[blockIdx.x*nThreads + threadIdx.x],blockIdx.x*nThreads + threadIdx.x,0,&states[blockIdx.x*nThreads + threadIdx.x]);
	int term = threadIdx.x + blockIdx.x*nThreads;
	int tPos = threadIdx.x*maxChain + blockIdx.x*maxChain*nThreads;
	for(w = 0; w<n_busca; w++)
	{
		do
		{
			for(i=0; i<inst->mAgents; i++)
			{
				res_aux[i]=0;
			}

			ejection->delta[term] = 0;
			ejection->op[term] = curand(&states[term])%2;
			__syncthreads();
			if(ejection->op[term] == 1)
			{
				ejection->sizeChain[term]=0;
				ejection->pos[0 + tPos] = curand(&states[term])%inst->nJobs;
				ejection->pos[1 + tPos] = curand(&states[term])%inst->mAgents;
				if((sol->resUsage[ejection->pos[1 + tPos] + blockIdx.x*inst->mAgents] + inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + ejection->pos[ 1 + tPos]] <= inst->capacity[ejection->pos[1 + tPos]])
						&&(tabuListshort[ejection->pos[0+tPos] + ejection->pos[1+tPos]*inst->nJobs + blockIdx.x*inst->nJobs*inst->mAgents]<iteration))
				{
					ejection->delta[term] = inst->cost[ejection->pos[0 + tPos]*inst->mAgents + ejection->pos[1 + tPos]] - inst->cost[ejection->pos[0 + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])];
					aux = 1;
				}
				//printf("Shift, with delta equals the %d\n",ejection->delta[term]);
			}
			else
			{

				aux = 1;
				ejection->sizeChain[term] = curand(&states[term])%(maxChain-1) + 2;
				size_aux = 2;
				d_best_aux = 1000000;
				ejection->pos[0 + tPos] = curand(&states[term])%inst->nJobs;
				ejection->delta[term] -= inst->cost[ejection->pos[0 + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0+tPos] + blockIdx.x*inst->nJobs])];
				ejection->pos[(ejection->sizeChain[term]-1) + tPos] = inst->nJobs-1;
				for(i=1; i<ejection->sizeChain[term]; i++)
				{
					ejection->pos[i + tPos] = curand(&states[term])%inst->nJobs;
					k = ejection -> pos[i + tPos];
					do
					{
						flag = 0;
						a = 0;

						for(j=0; j<i; j++)
						{
							if(ejection->pos[i + tPos]==ejection->pos[j + tPos])
							{
								flag = 1;
								break;
							}
						}

						if(((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs]) != ((int)sol->s[ejection->pos[(i-1) + tPos] + blockIdx.x*inst->nJobs]))
						{
							//						if(((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs])>4)
								//							printf("Agente: %d, Job: %d\n", ((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs]),ejection->pos[i + tPos]);
							if((flag!=1)&&(sol->resUsage[((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs]) + blockIdx.x*inst->mAgents] - inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs])] + inst->resourcesAgent[ejection->pos[(i-1) + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs])] <= inst->capacity[((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs])]))
							{
								if( tabuListshort[ ejection->pos[i-1 + tPos] + ((int)sol->s[ejection->pos[i + tPos] + blockIdx.x*inst->nJobs])*inst->nJobs + blockIdx.x*inst->nJobs*inst->mAgents]<iteration)
								{
									__syncthreads();
									if( i == ejection->sizeChain[term]-1)
									{
										if(((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs]) != ((int)sol->s[ejection->pos[(ejection->sizeChain[term]-1) + tPos] + blockIdx.x*inst->nJobs]))
										{
											if(tabuListshort[ ejection->pos[(ejection->sizeChain[term]-1) + tPos] + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])*inst->nJobs + blockIdx.x*inst->nJobs*inst->mAgents]<iteration )
											{
												res_aux[((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])] -= inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];
												res_aux[((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])] += inst->resourcesAgent[ejection->pos[(i-1) +tPos]*inst->mAgents + ((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];
												ejection->delta[term] += inst->cost[ejection->pos[(i-1) + tPos]*inst->mAgents+((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];//update delta
												ejection->delta[term] -= inst->cost[ejection->pos[i + tPos]*inst->mAgents+((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];

												if((sol->resUsage[((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs]) + blockIdx.x*inst->mAgents] - inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos]+blockIdx.x*inst->nJobs])] + inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])]<=inst->capacity[((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])])||(((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs]) == ((int)sol->s[ejection->pos[i + tPos]+blockIdx.x*inst->nJobs]))){
													delta_aux[i] = ejection->delta[term] + inst->cost[ejection->pos[i+tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])]; 
												}else{
													delta_aux[i] = 1000000;
												}
												if(delta_aux[i]<d_best_aux){
													d_best_aux = delta_aux[i];
													size_aux = i + 1;
												}
												a = 2;
												break; //if yes, next position is randomly selected
											}
										}
									}
									else
									{
										res_aux[((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])] -= inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];
										res_aux[((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])] += inst->resourcesAgent[ejection->pos[(i-1) +tPos]*inst->mAgents + ((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];
										ejection->delta[term] += inst->cost[ejection->pos[(i-1) + tPos]*inst->mAgents+((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];//update delta
										ejection->delta[term] -= inst->cost[ejection->pos[i + tPos]*inst->mAgents+((int)sol->s[ejection->pos[i +tPos] + blockIdx.x*inst->nJobs])];
										if((sol->resUsage[((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs]) + blockIdx.x*inst->mAgents] - inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos]+blockIdx.x*inst->nJobs])] + inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])]<=inst->capacity[((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])])||(((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs]) == ((int)sol->s[ejection->pos[i + tPos]+blockIdx.x*inst->nJobs]))){
											delta_aux[i] = ejection->delta[term] + inst->cost[ejection->pos[i+tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])];
										}else{
											delta_aux[i] = 1000000;
										}
										if(delta_aux[i]<d_best_aux){
											d_best_aux = delta_aux[i];
											size_aux = i + 1;
										}
										a = 2;
										break;
									}
								}
							}
						}
						ejection->pos[i + tPos] = (ejection->pos[i + tPos]+1)%(inst->nJobs);
					}
					while(ejection->pos[i + tPos]!=k);
					if(a==0)
					{
						if((i>1)&&(((int)sol->s[ejection->pos[0+tPos] + blockIdx.x*inst->nJobs]) != ((int)sol->s[ejection->pos[(i-1)+tPos]+blockIdx.x*inst->nJobs]))&&
								(tabuListshort[ ejection->pos[(i-1) + tPos] + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])*inst->nJobs + blockIdx.x*inst->nJobs*inst->mAgents]<iteration ))
						{
							ejection->sizeChain[term] = i;
						}
						else
						{
							aux = 0;

						}
						break;
					}
				}
				ejection->sizeChain[term] = size_aux;
				ejection->delta[term] = d_best_aux;
				//ejection->delta[term] += inst->cost[ejection->pos[(ejection->sizeChain[term]-1)+tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])];//update with last and first position
				res_aux[((int)sol->s[ejection->pos[0 +tPos] + blockIdx.x*inst->nJobs])] -= inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 +tPos] + blockIdx.x*inst->nJobs])];
				res_aux[((int)sol->s[ejection->pos[0 +tPos] + blockIdx.x*inst->nJobs])] += inst->resourcesAgent[ejection->pos[(ejection->sizeChain[term]-1) +tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 +tPos] + blockIdx.x*inst->nJobs])];
				if((sol->resUsage[((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs]) + blockIdx.x*inst->mAgents] - inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos]+blockIdx.x*inst->nJobs])] + inst->resourcesAgent[ejection->pos[(ejection->sizeChain[term]-1) +tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])]>inst->capacity[((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])])||(((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs]) == ((int)sol->s[ejection->pos[(ejection->sizeChain[term]-1) + tPos]+blockIdx.x*inst->nJobs])))
				{
					aux=0;
				}
				for(i=0; i<inst->mAgents; i++)
				{
					if(res_aux[i]>inst->capacity[i])
					{
						aux=0;
						break;
					}
				}
				//printf("Ejection Chain, with delta equals the %d\n",ejection->delta[term]);
			}

		}
		while(aux==0);
		
		if(w==0){
			op_best = ejection->op[term];
			delta_best = ejection->delta[term];
			for (i=0;i<maxChain;i++){
				pos_best[i] =  ejection->pos[i + tPos];
			}
			size_best = ejection->sizeChain[term];

		}else{
			if(ejection->delta[term]<delta_best){
				op_best = ejection->op[term];
				delta_best = ejection->delta[term];
				for (i=0;i<maxChain;i++){
					pos_best[i] =  ejection->pos[i + tPos];
				}
				size_best = ejection->sizeChain[term];
			}

		}
	}
	ejection->op[term] = op_best;
	ejection->delta[term] = delta_best;
	for (i=0;i<maxChain;i++){
		ejection->pos[i + tPos] = pos_best[i];
	}
	ejection->sizeChain[term] = size_best;
}






Solution* createGPUsolution(Solution* h_solution,TnJobs nJobs, TmAgents mAgents)
{
	size_t size_solution = sizeof(Solution)
                        		   + sizeof(TcostFinal)*nBlocks
                        		   + sizeof(Ts)*(nJobs*nBlocks) //vector s
                        		   + sizeof(TresUsage)*(mAgents*nBlocks); // vector resUsage



	Solution *d_sol;
	gpuMalloc((void**)&d_sol, size_solution);
	gpuMemset(d_sol,0,size_solution);
	h_solution->costFinal = (TcostFinal*)(d_sol+1);
	h_solution->s = (Ts*)(h_solution->costFinal + nBlocks);
	h_solution->resUsage = (TresUsage*)(h_solution->s + (nJobs*nBlocks));
	gpuMemcpy(d_sol, h_solution, size_solution, cudaMemcpyHostToDevice);
	return d_sol;
}

EjectionChain* createGPUejection(EjectionChain* h_ejection,TnJobs nJobs, TmAgents mAgents)
{
	size_t size_ejection = sizeof(EjectionChain)
                        		   + sizeof(Tpos)*(nBlocks*nThreads*maxChain)
                        		   + sizeof(Top)*(nBlocks*nThreads)
                        		   + sizeof(TSizeChain)*(nBlocks*nThreads)
                        		   + sizeof(Tdelta)*(nBlocks*nThreads);
	EjectionChain *d_ejection;
	gpuMalloc((void**)&d_ejection, size_ejection);
	gpuMemset(d_ejection,0,size_ejection);
	h_ejection->pos=(Tpos*)(d_ejection + 1);
	h_ejection->op = (Top*)(h_ejection->pos+ (nBlocks*nThreads*maxChain));
	h_ejection->sizeChain = (TSizeChain*)(h_ejection->op + (nBlocks*nThreads));
	h_ejection->delta = (Tdelta*)(h_ejection->sizeChain + (nBlocks*nThreads));
	gpuMemcpy(d_ejection, h_ejection, size_ejection, cudaMemcpyHostToDevice);
	return d_ejection;
}
