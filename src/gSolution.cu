#include "gSolution.cuh"

const int nBlocks =2;
const int nThreads = 100;
const int maxChain = 10;
//remove sizeTabu of parameters in gSolution.cuh and TS-GAP.cu

__global__ void TS_GAP(Instance *inst, Solution *sol, EjectionChain *ejection, int *tabuListshort, unsigned int *seed, curandState_t* states, int iteration, int n_busca){
	//variables of auxiliars
	int i, j,k,flag, aux, aux_2;
	//use for counting amount resources
	int res_aux[100];
	//Solution of block in memory shared
	__shared__ short int s_shared[1600];
//	s_shared = (short int *)malloc(sizeof(short int)*inst->nJobs);
	//Iterator of size search
	int s_search;
	int term = threadIdx.x + blockIdx.x*nThreads;
	int tPos = threadIdx.x*maxChain + blockIdx.x*maxChain*nThreads;

	//save best Delta of ejection chain
	int d_best_aux;
	//save best Delta of iteration
	int delta_best;
	//save best ejection chain
	short int pos_best[maxChain];
	short int size_best;
	short int op_best;
	//Delta for each step of the ejection chain 
	int delta_aux[maxChain];
	//save best size of ejection chain
	int size_aux;


	curand_init(seed[blockIdx.x*nThreads + threadIdx.x],blockIdx.x*nThreads + threadIdx.x,0,&states[blockIdx.x*nThreads + threadIdx.x]);
	aux = 0;
	aux_2 = inst->nJobs/nThreads;
	for(i = 0; i<= aux_2;i++){
		aux = inst->nJobs - i*nThreads; 
		if(threadIdx.x < aux){
			s_shared[threadIdx.x + i*nThreads] = ((int)sol->s[threadIdx.x + i*nThreads + blockIdx.x*inst->nJobs]);
			__syncthreads();
		}
		__syncthreads();
	}
	__syncthreads();
	
//	if((threadIdx.x == 0)&&(blockIdx.x==0)){
//		for(i=0;i<inst->nJobs;i++){
//			printf("job %d  agent %d\n", i,s_shared[i]);
//		}
		
//	}
	
	for(s_search = 0; s_search < n_busca; s_search++){
		do{
			aux = 0;
			for(i=0; i< inst->mAgents; i++){
				res_aux[i] = sol->resUsage[i + blockIdx.x*inst->mAgents]; 
			}
			ejection->delta[term] = 0;
			ejection->op[term] = curand(&states[term])%2;
			if(ejection->op[term] == 1){

				ejection->sizeChain[term] = 0;
				ejection->pos[0 + tPos] = curand(&states[term])%inst->nJobs;
				ejection->pos[1 + tPos] = curand(&states[term])%inst->mAgents;
				if((sol->resUsage[ejection->pos[1 + tPos] + blockIdx.x*inst->mAgents] + inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + ejection->pos[ 1 + tPos]] <= inst->capacity[ejection->pos[1 + tPos]])
						&&(tabuListshort[ejection->pos[0+tPos] + blockIdx.x*inst->nJobs]<=iteration))
				{
					ejection->delta[term] = inst->cost[ejection->pos[0 + tPos]*inst->mAgents + ejection->pos[1 + tPos]] - inst->cost[ejection->pos[0 + tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])];
					res_aux[ejection->pos[1 + tPos] + blockIdx.x*inst->mAgents] += inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + ejection->pos[ 1 + tPos]];
					res_aux[s_shared[ejection->pos[0 + tPos]]] -= inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]]];
					aux = 1;
				}
				d_best_aux = ejection->delta[term];
			}else{
				aux = 1;
				ejection->sizeChain[term] = maxChain;
				size_aux = 2;
				d_best_aux = 1000000;

				//choosing first job what not tabu list
				do{
					ejection->pos[0 + tPos] = curand(&states[term])%inst->nJobs;
				}while(tabuListshort[ejection->pos[0+tPos] + blockIdx.x*inst->nJobs]>iteration);

				ejection->delta[term] -= inst->cost[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0+tPos]]];
				ejection->pos[(ejection->sizeChain[term]-1) + tPos] = inst->nJobs-1;
				for(i=1; i<ejection->sizeChain[term]; i++){
					do{
						ejection->pos[i + tPos] = curand(&states[term])%inst->nJobs;
					}while(tabuListshort[ejection->pos[i+tPos] + blockIdx.x*inst->nJobs]>iteration);
					k = ejection -> pos[i + tPos];
					do{
						flag = 0;
						aux_2 = 0;
						//a = 0; //descobrir para que serve;
						for(j=0;  j < i ; j++){
							if(ejection->pos[i + tPos]==ejection->pos[j + tPos]){
								flag = 1;
								break;
							}
						}
						__syncthreads();
						if((flag!=1)&&(s_shared[ ejection->pos[i + tPos] ] != s_shared[ ejection->pos[(i-1) + tPos] ])){
							if( sol->resUsage[s_shared[ejection->pos[i + tPos]]+blockIdx.x*inst->mAgents] - inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i + tPos]]] + inst->resourcesAgent[ejection->pos[(i-1)+tPos]*inst->mAgents+s_shared[ejection->pos[i + tPos]]] <= inst->capacity[s_shared[ejection->pos[i + tPos] ]]){
								__syncthreads();
								if( i == ejection->sizeChain[term]-1){
									if((s_shared[ejection->pos[0 + tPos]])!=(s_shared[ejection->pos[(ejection->sizeChain[term]-1) + tPos]])){
										//res_aux[s_shared[ejection->pos[i +tPos]]] -= inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + s_shared[ejection->pos[i +tPos]]];
										res_aux[s_shared[ejection->pos[i +tPos]]] += inst->resourcesAgent[ejection->pos[(i-1) +tPos]*inst->mAgents + s_shared[ejection->pos[i +tPos]]] - inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + s_shared[ejection->pos[i +tPos]]];
										ejection->delta[term] += inst->cost[ejection->pos[(i-1) + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]] - inst->cost[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]]; //update delta
										//ejection->delta[term] -= inst->cost[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]];
										if(sol->resUsage[s_shared[ejection->pos[0 + tPos]] + blockIdx.x*inst->mAgents] - inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]]] + inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]]] <= inst->capacity[s_shared[ejection->pos[0 + tPos]]]){
											delta_aux[i] = ejection->delta[term] + inst->cost[ejection->pos[i+tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])]; 
										}else{
											delta_aux[i] = 1000000;
										}
										if(delta_aux[i]<d_best_aux){
											d_best_aux = delta_aux[i];
											size_aux = i + 1;
										}
										aux_2 = 2;
										break;
									}
								}else{
									//res_aux[s_shared[ejection->pos[i +tPos]]] -= inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + s_shared[ejection->pos[i +tPos]]];
									res_aux[s_shared[ejection->pos[i +tPos]]] += inst->resourcesAgent[ejection->pos[(i-1) +tPos]*inst->mAgents + s_shared[ejection->pos[i +tPos]]] - inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + s_shared[ejection->pos[i +tPos]]];
									ejection->delta[term] += inst->cost[ejection->pos[(i-1) + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]] - inst->cost[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]]; //update delta
									//ejection->delta[term] -= inst->cost[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]];
									if(sol->resUsage[s_shared[ejection->pos[0 + tPos]] + blockIdx.x*inst->mAgents] - inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]]] + inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]]] <= inst->capacity[s_shared[ejection->pos[0 + tPos]]]){
										delta_aux[i] = ejection->delta[term] + inst->cost[ejection->pos[i+tPos]*inst->mAgents + ((int)sol->s[ejection->pos[0 + tPos] + blockIdx.x*inst->nJobs])]; 
									}else{
										delta_aux[i] = 1000000;
									}
									if(delta_aux[i]<d_best_aux){
										d_best_aux = delta_aux[i];
										size_aux = i + 1;
									}
									aux_2 = 2;
									break;
								}
							}
						}	
						do{		
							ejection->pos[i + tPos] = (ejection->pos[i + tPos]+1)%(inst->nJobs);
						}while((tabuListshort[ejection->pos[i+tPos] + blockIdx.x*inst->nJobs]>iteration)&&(ejection->pos[i + tPos]!=k));
					}while(ejection->pos[i + tPos]!=k);
					if(aux_2==0)
                                        {
                                                if((i>1)&&(s_shared[ejection->pos[0+tPos]]) != (s_shared[ejection->pos[(i-1)+tPos]]))             
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
//                                res_aux[s_shared[ejection->pos[0 +tPos]]] -= inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0 +tPos]]];
                                res_aux[s_shared[ejection->pos[0 +tPos]]] += inst->resourcesAgent[ejection->pos[(ejection->sizeChain[term]-1) +tPos]*inst->mAgents + s_shared[ejection->pos[0 +tPos]]] -  inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0 +tPos]]];

                                if((sol->resUsage[s_shared[ejection->pos[0 + tPos]]] - inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]]] + inst->resourcesAgent[ejection->pos[(ejection->sizeChain[term]-1) +tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]]]>inst->capacity[s_shared[ejection->pos[0 + tPos]]])||(s_shared[ejection->pos[0 + tPos]] == s_shared[ejection->pos[(ejection->sizeChain[term]-1) + tPos]]))
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

			}

		}while((aux==0)||(d_best_aux== 1000000));
		if(s_search==0){
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
	__syncthreads();
//	if(threadIdx.x==0)
//		free(s_shared);
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

