#include  "guloso.h"

 


const int nThread =64;

//Create vector with priority of allocation jobs in Agents 
int* inicializeVector(Instance *inst, float p1, float p2)
{
    int *vOrdem;
    float *vParametro;
    int i, j;
    int aux1,aux2,iAux1,iAux2;
    vOrdem = (int*)malloc(sizeof(int)*(inst->nJobs *inst->mAgents));
    vParametro = (float*)malloc(sizeof(float)*(inst->nJobs *inst->mAgents));
    for(i=0; i<inst->nJobs ; i++)
    {
        for(j=0; j<inst->mAgents; j++)
        {
            vParametro[iReturn(i,j,inst->nJobs ,inst->mAgents)] = p1*inst->cost[iReturn(i,j,inst->nJobs ,inst->mAgents)] + p2*inst->resourcesAgent[iReturn(i,j,inst->nJobs ,inst->mAgents)];
            vOrdem[iReturn(i,j,inst->nJobs ,inst->mAgents)]=j;
        }
    }
    for(i=0; i<inst->nJobs ; i++)
    {
        for(j= inst->mAgents-1; j>=0; j--)
        {
            for(aux1= j-1; aux1>=0; aux1--)
            {
                if(vParametro[iReturn(i,j,inst->nJobs ,inst->mAgents)]<vParametro[iReturn(i,aux1,inst->nJobs ,inst->mAgents)])
                {
                    aux2 = vParametro[iReturn(i,j,inst->nJobs ,inst->mAgents)];
                    iAux2 = vOrdem[iReturn(i,j,inst->nJobs ,inst->mAgents)];
                    vParametro[iReturn(i,j,inst->nJobs ,inst->mAgents)]= vParametro[iReturn(i,aux1,inst->nJobs ,inst->mAgents)];
                    vOrdem[iReturn(i,j,inst->nJobs ,inst->mAgents)]=vOrdem[iReturn(i,aux1,inst->nJobs ,inst->mAgents)];;
                    vParametro[iReturn(i,aux1,inst->nJobs ,inst->mAgents)]=aux2;
                    vOrdem[iReturn(i,aux1,inst->nJobs ,inst->mAgents)]=iAux2;
                }
            }
        }
    }
    free(vParametro);
    return vOrdem;
}





int guloso(Instance *inst, Solution *sol, float p1, float p2, int block)
{
    int *vOrdem;
    int *allocated=(int*)malloc(sizeof(int)*inst->nJobs );
    int i,j,agent;
    int cont=0;
    memset(allocated,0,sizeof(int)*inst->nJobs );
    sol->costFinal[block]=0;
    vOrdem=inicializeVector(inst,p1,p2);
    for(i=0; i<inst->nJobs ; i++)
    {
        for(j=0; j<inst->mAgents; j++)
        {
        	//printf("res:%d\n",sol->resUsage[agent+block*inst->mAgents]);
            agent=vOrdem[iReturn(i,j,inst->nJobs ,inst->mAgents)];
            /*printf("n: %d \n",iReturn(i,j,inst->nJobs ,inst->mAgents));*/
            if((allocated[i]==0)&&(inst->resourcesAgent[iReturn(i,agent,inst->nJobs ,inst->mAgents)]+sol->resUsage[agent + block*inst->mAgents]<=inst->capacity[agent]))
            {
                allocated[i]=1;
                sol->s[i + block*inst->nJobs] = ((char)agent);
                sol->costFinal[block]+=inst->cost[iReturn(i,agent,inst->nJobs ,inst->mAgents)];
                sol->resUsage[agent + block*inst->mAgents]+=inst->resourcesAgent[iReturn(i,agent,inst->nJobs ,inst->mAgents)];
                cont++;
            }
        }
    }

    if(cont!=inst->nJobs ){
        sol->costFinal[block]=0;
        return 0;
    }
    return 1;
    free(allocated);
    free(vOrdem);

}

