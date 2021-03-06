#include "gSolution.cuh"

#define nBlocks 10 //const int nBlocks =10;
#define nThreads 64 //const int nThreads = 64;
#define maxChain 10 //const int maxChain = 10; // tamanho maximo de uma cadeia
//remove sizeTabu of parameters in gSolution.cuh and TS-GAP.cu

__global__ void TS_GAP(Instance *inst, Solution *sol, EjectionChain *ejection, int *tabuListshort, unsigned int *seed, curandState_t* states, int iteration, int n_busca){
        // n_busca eh o numero de CADEIAS que uma thread vai explorar

	//variables of auxiliars
	register int i, j,k,flag, aux, aux_2;
        register short int nJobs = inst->nJobs;
        register short int mAgents = inst->mAgents;

	//use for counting amount resources
	register unsigned short int res_aux[800];
	//Solution of block in memory shared
	//__shared__ short int s_shared[9000];
	__shared__ unsigned char s_shared[16000];
	//s_shared = (short int *)malloc(sizeof(short int)*inst->nJobs);
	//Iterator of size search
	int s_search;
	#define term  threadIdx.x + blockIdx.x*nThreads
	#define tPos  threadIdx.x*maxChain + blockIdx.x*maxChain*nThreads

	//save best Delta of ejection chain
	int d_best_aux;
	//save best Delta of iteration
	int delta_best;
	//save best ejection chain
	register short int pos_best[maxChain];
	register short int size_best;
	register short int op_best;
	//Delta for each step of the ejection chain 
	register int delta_aux[maxChain]; // substituir por melhor_delta, melhor_pos
	//save best size of ejection chain
	int size_aux;


	curand_init(seed[blockIdx.x*nThreads + threadIdx.x],blockIdx.x*nThreads + threadIdx.x,0,&states[blockIdx.x*nThreads + threadIdx.x]);
	aux = 0;
	//aux_2 = nJobs/nThreads;
	for(i = 0; i<= nJobs/nThreads; ++i){ // <= aux_2
		//aux = nJobs - i*nThreads; 
		if(threadIdx.x <  nJobs - i*nThreads){ // < aux
			s_shared[threadIdx.x + i*nThreads] = sol->s[threadIdx.x + i*nThreads + blockIdx.x*inst->nJobs];
			//__syncthreads();
		}
		//__syncthreads();
	}
	__syncthreads();
	
//	if((threadIdx.x == 0)&&(blockIdx.x==0)){
//		for(i=0;i<inst->nJobs;i++){
//			printf("job %d  agent %d\n", i,s_shared[i]);
//		}		
//	}

        // numero de cadeias que cada thread ira fazer
	for(s_search = 0; s_search < n_busca; ++s_search){
		do{
			aux = 0; // flag de movimento infactivel
			for(i=0; i < mAgents; ++i){
				res_aux[i] = sol->resUsage[i + blockIdx.x*mAgents]; 
			}
			ejection->delta[term] = 0;
			ejection->op[term] = curand(&states[term])%2;

                        // op=0, ejection chain
                        // op=1, troca
			if(ejection->op[term] == 1){ 
                                // TROCA
				ejection->sizeChain[term] = 0;
				register short int swap_job = curand(&states[term])%nJobs;   // ex: 8500 (agente 30)
				register unsigned char swap_agent = curand(&states[term])%mAgents; // ex: 90
                                ejection->pos[0 + tPos] = swap_job; 
                                ejection->pos[1 + tPos] = swap_agent; 

                                // USA POSIÇÕES 0 e 1 para troca
                                // consigo tirar a job 8500 do agente 30 e colocar no agente 90?

                                // ejection->pos:  th0 bl0 |0,1,2, ..., maxChain|; th1 bl0 |maxChain+1, maxChain+2 ... 2*maxChain| ... 
				if((sol->resUsage[swap_agent + blockIdx.x*mAgents] + inst->resourcesAgent[swap_job*mAgents + swap_agent] <= inst->capacity[swap_agent])
						&&(tabuListshort[swap_job + blockIdx.x*nJobs]<=iteration))
				{
					ejection->delta[term] = inst->cost[swap_job*mAgents + swap_agent] - inst->cost[swap_job*mAgents + s_shared[swap_job]];

					res_aux[swap_agent] += inst->resourcesAgent[swap_job*mAgents + swap_agent];
					res_aux[s_shared[swap_job]] -= inst->resourcesAgent[swap_job*mAgents + s_shared[swap_job]];
					aux = 1; // movimento eh factivel
				}
				d_best_aux = ejection->delta[term];
			}else{
                                // op == 0, EJECTION CHAIN
				aux = 1; // flag de movimento factivel
				ejection->sizeChain[term] = maxChain;
				size_aux = 2;
				d_best_aux = 1000000;

				//choosing first job what not tabu list
                                register short int first_job = 0; 
				do{
					first_job = curand(&states[term])%nJobs; // sorteou job 8500
                                        // job escolhido. Sera que eh tabu?
				}while(tabuListshort[first_job + blockIdx.x*nJobs] > iteration);
				ejection->pos[0 + tPos] = first_job;

				ejection->delta[term] -= inst->cost[first_job*mAgents + s_shared[first_job]]; // s_shared[8500] => agente 30
				ejection->pos[(maxChain-1) + tPos] = nJobs-1; // ????
				for(i=1; i<maxChain; ++i){
					do{
						k = curand(&states[term])%inst->nJobs;
                                                // job escolhido. Sera que ele eh tabu?
					}while(tabuListshort[k + blockIdx.x*nJobs]>iteration);
					ejection->pos[i + tPos] = k;
                                        
					do{
						flag = 0; // indica que nao ha repeticao
						aux_2 = 0; // indica que o ultimo passo da ejection chain eh invalido
						//a = 0; //descobrir para que serve;
						for(j=0;  j < i ; j++){
							if(ejection->pos[i + tPos]==ejection->pos[j + tPos]){
								flag = 1; // indica que ha repeticao
								break;
							}
						}
						/////////__syncthreads();

                                                // jobs diferentes (flag==0), mas testa agora se os AGENTES (s_shared) sao diferentes...
						if((flag==0) && (s_shared[ ejection->pos[i + tPos] ] != s_shared[ ejection->pos[(i-1) + tPos] ])){
							if( sol->resUsage[s_shared[ejection->pos[i + tPos]]+blockIdx.x*mAgents] - inst->resourcesAgent[ejection->pos[i + tPos]*mAgents+s_shared[ejection->pos[i + tPos]]] + inst->resourcesAgent[ejection->pos[(i-1)+tPos]*mAgents+s_shared[ejection->pos[i + tPos]]] <= inst->capacity[s_shared[ejection->pos[i + tPos] ]]){
								//__syncthreads();
								if(i == maxChain-1){  //if( i == ejection->sizeChain[term]-1){
									if(s_shared[first_job]!=s_shared[ejection->pos[(maxChain-1) + tPos]]){
										//res_aux[s_shared[ejection->pos[i +tPos]]] -= inst->resourcesAgent[ejection->pos[i + tPos]*inst->mAgents + s_shared[ejection->pos[i +tPos]]];
										res_aux[s_shared[ejection->pos[i +tPos]]] += inst->resourcesAgent[ejection->pos[(i-1) +tPos]*mAgents + s_shared[ejection->pos[i +tPos]]] - inst->resourcesAgent[ejection->pos[i + tPos]*mAgents + s_shared[ejection->pos[i +tPos]]];
										ejection->delta[term] += inst->cost[ejection->pos[(i-1) + tPos]*mAgents+s_shared[ejection->pos[i +tPos]]] - inst->cost[ejection->pos[i + tPos]*mAgents+s_shared[ejection->pos[i +tPos]]]; //update delta
										//ejection->delta[term] -= inst->cost[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]];
										if(sol->resUsage[s_shared[first_job] + blockIdx.x*mAgents] - inst->resourcesAgent[first_job*mAgents + s_shared[first_job]] + inst->resourcesAgent[ejection->pos[i + tPos]*mAgents + s_shared[ejection->pos[first_job]]] <= inst->capacity[s_shared[first_job]]){
											delta_aux[i] = ejection->delta[term] + inst->cost[ejection->pos[i+tPos]*mAgents + s_shared[first_job] ]; 
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
									res_aux[s_shared[ejection->pos[i +tPos]]] += inst->resourcesAgent[ejection->pos[(i-1) +tPos]*mAgents + s_shared[ejection->pos[i +tPos]]] - inst->resourcesAgent[ejection->pos[i + tPos]*mAgents + s_shared[ejection->pos[i +tPos]]];
									ejection->delta[term] += inst->cost[ejection->pos[(i-1) + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]] - inst->cost[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]]; //update delta
									//ejection->delta[term] -= inst->cost[ejection->pos[i + tPos]*inst->mAgents+s_shared[ejection->pos[i +tPos]]];
									if(sol->resUsage[s_shared[first_job] + blockIdx.x*mAgents] - inst->resourcesAgent[first_job*mAgents + s_shared[first_job]] + inst->resourcesAgent[ejection->pos[i + tPos]*mAgents + s_shared[first_job]] <= inst->capacity[s_shared[first_job]]){
										delta_aux[i] = ejection->delta[term] + inst->cost[ejection->pos[i+tPos]*inst->mAgents + s_shared[ejection->pos[0 + tPos]] ]; 
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

                                                // garante que o elemento sorteado para posicao i sera diferente dos anteriores (busca sequencial) 
						do{		
							ejection->pos[i + tPos] = (ejection->pos[i + tPos]+1)%nJobs;
						} while((tabuListshort[ejection->pos[i+tPos] + blockIdx.x*nJobs] > iteration)&&(ejection->pos[i + tPos]!=k));

                                                // se 'der a volta', ele ja sai do proximo while tambem e desiste de aumentar a ejection chain
					} while(ejection->pos[i + tPos]!=k);

                                        // se o ultimo passo da ejection chain for invalido, entre no if
					if(aux_2==0)
                                        {
                                                // sera que primeiro e ultimo estao no mesmo agente? (isto eh possivel??)
                                                if( (i>1) && (s_shared[first_job] != s_shared[ejection->pos[(i-1)+tPos]]) )             
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
                                //res_aux[s_shared[ejection->pos[0 +tPos]]] -= inst->resourcesAgent[ejection->pos[0 + tPos]*inst->mAgents + s_shared[ejection->pos[0 +tPos]]];
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

