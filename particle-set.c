#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define NSIZE 50000     // path node table size
#define LSIZE 50000     // path line table size
#define PSIZE 9999999   // raindrop particle table size
#define MCELL 200       // the cell size of the flow map

#define MEANV 0.0675
#define DENSITY 10

typedef struct{
    int nPNID;
    
    double fX;
    double fY;
    double fZ;
    
    double fin;     // number of in lines
    double fout;    // number of out lines
    
    int ntype;      // if a source node
    int nPPRID;     // assign to a precipitation point
    
    int nPLID_2;    // from nPNID to path line nPLID_2
    int nPNID_2;    // from nPNID to path node nPNID_2
}PathNode;

typedef struct{
    int nPLID;
    
    int nPNID_1;
    int nPNID_2;
    int nPLID_2;    // connected path line
    
    int ntype;      // if from a source node
    
    double fLength;
    double fVelocity;
    double fTime;
}PathLine;

typedef struct{
    PathNode *pNode_List; // path node list
    PathLine *pLine_List; // path line list
    
    int *pNid_List; // node id list
    int *pSNid_List; // source node id list
    int *pBNid_List; // basin node id list
    
    int *pLid_List; // line id list
    
    int nNumNodes;
    int nNumSources;
    int nNumBasins;
    
    int nNumLines;
    
    double fLeft;
    double fRight;
    double fTop;
    double fBottom; // the rect range of the drainge area
}DrainageNet;

typedef struct{
    int nPID;
    
    int nPLID_on; // the path line inside
    int nstatus; // whether moving or stopped
    double fPos; // relative position on the path line
    
    double fX;
    double fY; // the geographic coordinates
}Particle;

/* MPI struct for loading data */
void new_mpi_struct(MPI_Datatype* mpi_particle) {
    int blockcounts[2] = { 3, 3};						// two data blocks
    MPI_Aint offsets[2] = { 0, 3 * sizeof(int) };		// completion with 4 bytes
    MPI_Datatype oldtypes[2] = { MPI_INT, MPI_DOUBLE };	// block data type
    
    MPI_Type_struct(2, blockcounts, offsets, oldtypes, mpi_particle);
    MPI_Type_commit(mpi_particle);
}

typedef struct{
    int *pSN_List; // selected source node id list
    int *pRD_List; // the number of the raindrops into each source node
    
    int nNumSSN; // number of the selected source nodes
}Rain;

typedef struct{
    int nPPRID;
    
    double *pfPPR; // the precipitation rate mm/hr
    int nTimeFrame; // time frames
}PPRate; //precipitation rate points

typedef struct{
    PPRate *pPPR_List;
    int nNumPoints; // the number of grid cells
    int nNumFrames; // the number of time frames
}PPMap; //dynamic precipitation map

typedef struct{
    Particle *pRP_List; // raining particle list
    int nNumParticles; // the number of particles
}Runoff;

typedef struct{
    double fSize; // length and width of each cell
    int nRows; // number of cells (north-to-sourth)
    int nColumns; // number of cells (west-to-east)
    int *pFlowMap; // 2D matrix for number of particles in each cell
}FlowMap;

// merging strings abc and ABC to abcABC
char* join_string(char *s1, char *s2) {
    char *result = malloc(strlen(s1)+strlen(s2)+1);
    if (result == NULL) exit (1);
    
    strcpy(result, s1);
    strcat(result, s2);
    
    return result;
}

PPMap * loading_rainfall(char *pFilePath) {
    FILE *fpPPR;
    
    int nNum=0, nID, nPts, nFrames, i;
    double ftmp;
    
    PPMap *pPPMap = (PPMap*) malloc(sizeof(PPMap));
    
    char* pPPRPath = join_string(pFilePath, "/input/rainfall.in");
    
    if((fpPPR = fopen(pPPRPath, "r"))==NULL) {
        printf("precipitation file open error\n");
        free(pPPRPath);
        pPPRPath = NULL;
        return pPPMap;
    } else {
        fscanf(fpPPR,"%d\t%d\n", &nPts, &nFrames);
        
        pPPMap->nNumPoints = nPts;
        pPPMap->nNumFrames = nFrames;
        pPPMap->pPPR_List = (PPRate*) malloc(sizeof(PPRate)*nPts);
        
        while (!feof(fpPPR)) {
            
            fscanf(fpPPR,"%d\t", &nID);
            pPPMap->pPPR_List[nNum].nPPRID = nID;
            
            pPPMap->pPPR_List[nNum].pfPPR = (double*) malloc(sizeof(double)*pPPMap->nNumFrames);
            for (i=0; i<pPPMap->nNumFrames; i++) {
                fscanf(fpPPR,"%lf\t", &ftmp);
                pPPMap->pPPR_List[nNum].pfPPR[i] = ftmp;
            }
            fscanf(fpPPR,"\n");
            
            nNum++;
        }
        
        pPPMap->nNumPoints = nNum;
        //printf("rainfall points:%d\nrainfall %d hours\n", nNum, pPPMap->nNumFrames);
        free(pPPRPath);
        pPPRPath = NULL;
        
        fclose(fpPPR);
    }
    
    return pPPMap;
}

int GetIndexByID(int nID, PPMap * pPPMap) {
    int i;
    for (i=0; i<pPPMap->nNumPoints; i++) {
        if (pPPMap->pPPR_List[i].nPPRID == nID) {
            int nppr = i;
            return nppr;
        }
    }
    return -1;
}

// calculating length of flow path line
double get_length(PathNode node_1, PathNode node_2) {
    return sqrtf((node_1.fX-node_2.fX)*(node_1.fX-node_2.fX) + (node_1.fY-node_2.fY)*(node_1.fY-node_2.fY) + (node_1.fZ-node_2.fZ)*(node_1.fZ-node_2.fZ));
}

// relating path nodes and lines with flow direction
int topology(DrainageNet *net) {
    
    PathNode *node = (PathNode*) malloc(sizeof(PathNode));
    PathNode *node_to = (PathNode*) malloc(sizeof(PathNode));
    PathLine *line_on = (PathLine*) malloc(sizeof(PathLine));
    
    int nsn = 0, nbn = 0, npn = 0, njn = 0;
    int id_node, id_line, i;
    
    // relating nodes
    for (i=0; i<net->nNumNodes; i++) {
        *node = net->pNode_List[net->pNid_List[i]];
        
        id_node = -1;
        id_line = -1;
        
        if (node->fin != 0.0 && node->fout == 0.0) {
            net->pBNid_List[nbn++] = net->pNid_List[i];
            net->pNode_List[node->nPNID].ntype = 2; // outlet
        } else if (node->fin == 0.0 && node->fout != 0.0) {
            net->pSNid_List[nsn++] = net->pNid_List[i];
            net->pNode_List[node->nPNID].ntype = 0; // source
        }  else if (node->fin != 0.0 && node->fout != 0.0) {
            npn ++;
            net->pNode_List[node->nPNID].ntype = 1; // passing
        } else {
            njn++;  // junction
        }
    }
    
    net->nNumSources = nsn;
    net->nNumBasins = nbn;
    
    // relating lines
    int nsl = 0, npl = 0, nbl = 0, nsbl = 0;
    for (i=0; i<net->nNumLines; i++) {
        *line_on = net->pLine_List[net->pLid_List[i]];
        
        *node = net->pNode_List[line_on->nPNID_1];
        *node_to = net->pNode_List[line_on->nPNID_2];
        
        net->pLine_List[line_on->nPLID].nPLID_2 = node_to->nPLID_2;
        //net->pLine_List[line_on.nPLID].fTime = get_time(node, node_to);
        
        if (node->ntype == 0 && node_to->ntype == 2) {
            nsbl += 1;
            net->pLine_List[line_on->nPLID].ntype = 3; // from source to basin
        } else if (node->ntype == 0 && node_to->ntype != 2) {
            nsl += 1;
            net->pLine_List[line_on->nPLID].ntype = 0; // from source
        } else if (node->ntype != 0 && node_to->ntype == 2) {
            nbl += 1;
            net->pLine_List[line_on->nPLID].ntype = 2; // to basin
            //printf("%d\n", line_on->nPLID);
        } else {
            npl += 1;
            net->pLine_List[line_on->nPLID].ntype = 1; // passing lines
        }
        
        //printf("%d,%f\n", line_on.nPLID, net->pLine_List[line_on.nPLID].fTime);
    }
    
    free(node);node=NULL;
    free(node_to);node_to=NULL;
    free(line_on);line_on=NULL;
    
    //printf("lines:%d\nline from source:%d\npassing line:%d\nline to basin:%d\nline from source to basin:%d\n\n", net->nNumLines, nsl, npl, nbl, nsbl);
    
    //printf("nodes:%d\nsources:%d\njunctions:%d\npassing:%d\nbasins:%d\n\n", net->nNumNodes, net->nNumSources, njn, npn, net->nNumBasins);
    
    if (net->nNumSources < 1 || net->nNumBasins < 1) {
        //printf("\nFailed to constructing drainge topology!\n\n");
        return 0;
    } else {
        //printf("\ndrainge topology construction completed!\n\n");
        //printRoute(net);
        return 1;
    }
}


// relating the flow in and out degrees of path nodes
int geometry(DrainageNet *net) {
    PathNode* node;
    
    // print the degrees of path nodes
    int *nDegrees = (int*)malloc(sizeof(int*)*225);
    
    int id_1, id_2, i, j;
    
    for (i=0; i<225; i++) {
        nDegrees[i] = 0;
    }
    
    for (i=0; i<net->nNumNodes; i++) {
        node = &net->pNode_List[net->pNid_List[i]];
        id_1 = -1;
        id_2 = -1;
        
        nDegrees[(int)(net->pNode_List[node->nPNID].fin*2)*15+(int)(net->pNode_List[node->nPNID].fout*2)] += 1;
    }
    
    //printf("\n********************\n");
    for (i=0;i<15;i+=2) {
        for (j=0;j<15;j+=2) {
            //printf("%d\t",nDegrees[i*15+j]);
        }
        //printf("\n");
    }
    //printf("********************\n\n");
    
    int nIfCorrect = 1;
    
    for (i=0; i<15; i++) {
        for (j=4; j<15; j++) {
            if (nDegrees[i*15+j] > 0) {
                nIfCorrect = 0;
            }
        }
    }
    
    free(nDegrees);
    nDegrees =NULL;
    
    if (!nIfCorrect) {
        //printf("\nCorrecting drainge geometry...\n\n");
        geometry(net);
        return 0;
    } else {
        //printf("\nGeometry correction completed!\n\n");
        return 1;
    }
}


int creating_net(char *pFilePath, DrainageNet *net) {
    
    FILE *fpN, *fpTopo;
    
    char str[240];
    int ntmp, nidn, nidl, ntmp2;
    double fx, fy, fz, ftemp;//, ftemp2;
    int nn = 0, nl = 0; // for numbers of nodes and lines
    
    PathNode *node_1 = (PathNode*) malloc(sizeof(PathNode));
    PathNode *node_2 = (PathNode*) malloc(sizeof(PathNode));
    
    char* pNPath = join_string(pFilePath, "/input/pathnode.in");
    if((fpN = fopen(pNPath, "r"))==NULL) {
        //printf("path node file open error\n");
        free(pNPath);pNPath=NULL;
        return -1;
    } else {
        fscanf(fpN,"%s\n",str);
        //printf("%s\n",str);
        
        while (!feof(fpN)) {
            fscanf(fpN,"%d, %lf, %lf, %lf\n", &nidn, &fx, &fy, &fz);
            
            net->pNode_List[nidn].nPNID = nidn;
            net->pNode_List[nidn].fX = fx;
            net->pNode_List[nidn].fY = fy;
            net->pNode_List[nidn].fZ = fz;
            
            net->pNode_List[nidn].nPPRID = nidn;
            
            net->pNid_List[nn++] = nidn;
        }
        
        net->nNumNodes = nn;
        //printf("flow path nodes:%d\n\n", nn);
        fclose(fpN);
        free(pNPath);
        pNPath =NULL;
    }
    
    
    // topology pathnode and pathline
    // char* pTopoPath = join_string(pFilePath, "/topology_sdz.txt");
    char* pTopoPath = join_string(pFilePath, "/input/pathline.in");
    if((fpTopo = fopen(pTopoPath, "r"))==NULL) {
        printf("topology file open error\n");
        free(pTopoPath);
        pTopoPath=NULL;
        return -1;
    } else {
        fscanf(fpTopo,"%s",str);
        //printf("%s\n",str);
        while (!feof(fpTopo)) {
            
            fscanf(fpTopo,"%d, %d, %d, %lf\n", &nidl, &ntmp, &ntmp2, &ftemp);
            
            net->pLine_List[nidl].nPNID_1 = ntmp;
            net->pLine_List[nidl].nPNID_2 = ntmp2;
            
            if (ntmp2 != -1) {
                *node_1 = net->pNode_List[ntmp];
                *node_2 = net->pNode_List[ntmp2];
                
                // calculating the flow time on path
                net->pLine_List[nidl].fLength = get_length(*node_1, *node_2);
                
                net->pLine_List[nidl].fVelocity = 60*MEANV*ftemp/0.06; // from m/s to m/min
                
                net->pLine_List[nidl].fTime = get_length(*node_1, *node_2)/net->pLine_List[nidl].fVelocity;
                
                net->pNode_List[ntmp].nPLID_2 = nidl;
                net->pNode_List[ntmp].nPNID_2 = ntmp2;
                
                net->pNode_List[ntmp].fout += 1.0;
                net->pNode_List[ntmp2].fin += 1.0;
                
            } else {
                printf("wrong topology: isolate node: %d\n", ntmp);
            }
            
            net->pLid_List[nl++] = nidl;
        }
        net->nNumLines = nl;
        //printf("flow path lines:%d\n\n", nl);
        
        fclose(fpTopo);
        free(pTopoPath);
        pTopoPath =NULL;
    }
    
    free(node_1);node_1=NULL;
    free(node_2);node_2=NULL;
    
    if (net->nNumNodes < 1 || net->nNumLines < 1) {
        //printf("\nFailed to initializing drainge network!\n\n");
        return 0;
    } else {
        //printf("\ndrainge network initialization completed!\n\n");
        return 1;
    }
}

// preparing drainage network and time network
DrainageNet *initializing(int nNSize, int nLSize) {
    
    DrainageNet *net = (DrainageNet*) malloc(sizeof(DrainageNet));
    
    net->pNode_List = (PathNode*) malloc(sizeof(PathNode)*nNSize);
    net->pLine_List = (PathLine*) malloc(sizeof(PathLine)*nLSize);
    
    net->pNid_List = (int*) malloc(sizeof(int)*nNSize);
    net->pSNid_List = (int*) malloc(sizeof(int)*nNSize);
    net->pBNid_List = (int*) malloc(sizeof(int)*nNSize);
    
    net->pLid_List = (int*) malloc(sizeof(int)*nLSize);
    
    net->nNumNodes = 0;
    net->nNumLines = 0;
    net->nNumSources = 0;
    net->nNumBasins = 0;
    
    net->fLeft = 9999999.9;
    net->fRight = -9999999.9;
    net->fTop = -9999999.9;
    net->fBottom = 9999999.9;
    
    int i;
    
    for(i=0; i<nNSize; i++) {
        net->pNode_List[i].nPNID = i;
        net->pNode_List[i].fX = 0;
        net->pNode_List[i].fY = 0;
        net->pNode_List[i].fZ = 0;
        
        net->pNode_List[i].fin = 0;
        net->pNode_List[i].fout = 0;
        
        net->pNode_List[i].ntype = -1;
        net->pNode_List[i].nPPRID = -1;
        
        net->pNode_List[i].nPLID_2 = -1;
        net->pNode_List[i].nPNID_2 = -1;
        
        net->pNid_List[i] = -1;
        net->pSNid_List[i] = -1;
        net->pBNid_List[i] = -1;
    }
    
    for(i=0; i<nLSize; i++) {
        
        net->pLine_List[i].nPLID = i;
        
        net->pLine_List[i].nPNID_1 = -1;
        net->pLine_List[i].nPNID_2 = -1;
        net->pLine_List[i].nPLID_2 = -1;
        
        net->pLine_List[i].ntype = -1;
        net->pLine_List[i].fTime = 0;
        
        net->pLine_List[i].fLength = 0.0;
        net->pLine_List[i].fVelocity = 0.0;
    }
    
    return net;
}

// updating postition for a single particle
int moving(DrainageNet *net, double fTime, Particle *particle) {
    double fRemain = fTime;
    int nEnd = 0;
    PathLine *line = (PathLine*) malloc(sizeof(PathLine));
    PathNode *node_1 = (PathNode*) malloc(sizeof(PathNode));
    PathNode *node_2 = (PathNode*) malloc(sizeof(PathNode));
    
    *line = net->pLine_List[particle->nPLID_on];
    
    // to the next line
    while (fRemain >= line->fTime*(1-particle->fPos)) {
        if (line->nPLID_2 != -1) {
            
            fRemain = fRemain - line->fTime*(1-particle->fPos);
            
            particle->fPos = 0.0;
            particle->nPLID_on = line->nPLID_2;
            *line = net->pLine_List[particle->nPLID_on];
            
            // reach the basin outlet
            if (line->ntype == 2) {
                nEnd = 1;
            }
        } else {
            fRemain = 0.0; // reach basin
            
            particle->fPos = 1.0;
            particle->nstatus = 1;
            
            //nEnd = 0;
            break;
        }
        
    }
    
    // to the current line
    particle->fPos = particle->fPos + fRemain/line->fTime;
    
    // updating geographic coordinates
    *node_1 = net->pNode_List[line->nPNID_1];
    *node_2 = net->pNode_List[line->nPNID_2];
    particle->fX = node_1->fX + particle->fPos*(node_2->fX - node_1->fX);// - net->fLeft;
    particle->fY = node_1->fY + particle->fPos*(node_2->fY - node_1->fY);// - net->fBottom;
    
    free(line);line=NULL;
    free(node_1);node_1=NULL;
    free(node_2);node_2=NULL;
    
    return nEnd;
}

int main(int argc, char *argv[]) {
    /* setup the communication world */
    int my_rank, comm_sz, proc_len;					// parallel processors
    char proc_name[MPI_MAX_PROCESSOR_NAME];
    
    MPI_Init(&argc, &argv);							// initialize
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);		// num of processors
    MPI_Get_processor_name(proc_name, &proc_len);	// processor names
    
    //printf("%d\t%s\n",comm_sz, proc_name);
    
    if (comm_sz < 1) {
        fprintf(stderr,"No enough processors.\n");
        exit(-1);
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);		// processor ranks
    
    /*** for all processors***/
    clock_t start, finish, start_export, finish_export;
    double fDuration = 0, fDuration_export = 0;
    
    // drainge network
    FILE *fpM;
    int nx, ny;
    char strFilePath[120] = "/home/ubuntu/aws";
    char* pMPath = join_string(strFilePath, "/output/flowmap.out");
    char* pOPath = join_string(strFilePath, "/output/outlet.out");
    char* pTPath = join_string(strFilePath, "/output/timecost.out");
    
    // checking
    DrainageNet *net = initializing(NSIZE, LSIZE);
    creating_net(strFilePath, net);
    geometry(net);
    topology(net);
    
    // loading data (precipitation, infiltration, evapotranspiration)
    PPMap *pPPMap = loading_rainfall(strFilePath);
    
    // raining and runing
    int nST = 1;    // starting time of raining
    int nIT = 1;    // time interval (1 hour = 60 mins)
    int nET = 144;  // ending time of raining
    int nMP = 0;    // number of moving/outleting particles
    
    Runoff *runoff = (Runoff*) malloc(sizeof(Runoff));
    runoff->pRP_List = (Particle*) malloc(sizeof(Particle)*PSIZE);
    runoff->nNumParticles = 0;
    
    FlowMap *map = (FlowMap*) malloc(sizeof(FlowMap));
    net->fTop = 4232200;
    net->fBottom = 4217200;
    net->fLeft = 2524800;
    net->fRight = 2535200;
    map->fSize = MCELL;
    map->nRows = 1 + (net->fTop - net->fBottom)/map->fSize;
    map->nColumns = 1 + (net->fRight - net->fLeft)/map->fSize;
    map->pFlowMap =  (int*) malloc(sizeof(int)*map->nRows*map->nColumns);
    
    int ndest, nppr, nT, i, j, nPPRate;
    
    start = clock();
    for(nT=nST; nT<=nET; ) {
        nPPRate = 0; // rainfall particles (mm)
        nMP = 0;
        
        // 1. updating existing particles
        
        for(i=0; i<runoff->nNumParticles; i++) {
            if (runoff->pRP_List[i].nstatus == 0) { // still moving
                
                Particle *particle = (Particle*) malloc(sizeof(Particle));
                *particle = runoff->pRP_List[i];
                
                if (moving(net, 1*60.0, particle)) {
                    nMP++;
                }
                
                runoff->pRP_List[i] = *particle;
                free(particle);
                particle = NULL;
            }
        }
        
        // 2. adding new raindrop particles
        for(i=0; i<net->nNumSources; i++) {
            PathNode *node = (PathNode*) malloc(sizeof(PathNode));
            *node = net->pNode_List[net->pSNid_List[i]];
            
            ndest = (i+1)%(comm_sz);
            if (my_rank == ndest) {
                
                if (node->ntype == 0 && node->nPPRID != -1) { // the source nodes
                    
                    nppr = GetIndexByID(node->nPPRID, pPPMap);
                    
                    if (nppr == -1) {
                        nPPRate = 0;
                    } else {
                        nPPRate = (int) DENSITY*pPPMap->pPPR_List[nppr].pfPPR[nT-1];
                        
                        for(j=0; j<nPPRate; j++) {
                            Particle *particle = (Particle*) malloc(sizeof(Particle));
                            particle->nPID = runoff->nNumParticles;
                            particle->nPLID_on = node->nPLID_2;
                            particle->fPos = 0.0;
                            particle->nstatus = 0;
                            
                            if (moving(net, 1*60.0-(double)j*60.0/(double)nPPRate, particle)) {
                                nMP++;
                            }
                            
                            runoff->pRP_List[runoff->nNumParticles++] = *particle;
                            free(particle);
                            particle = NULL;
                        }
                    }
                }
            }
            
            free(node);
            node = NULL;
        }
        
        // export flow maps and outlet discharges
        start_export = clock();
        if (comm_sz == 1) {
            
            Particle *particle = (Particle*) malloc(sizeof(Particle));
            
            for(i=0;i<map->nRows;i++) {
                for(j=0;j<map->nColumns;j++) {
                    map->pFlowMap[i*map->nColumns+j] = 0;
                }
            }
            
            for (i=0; i<runoff->nNumParticles; i++) {
                *particle = runoff->pRP_List[i];
                
                //printf("%f %f %f\n", particle->fX, particle->fY, net->fLeft);
                if ((particle->fX - net->fLeft)/map->fSize < 0) {
                    nx = 0;
                } else if ((particle->fX - net->fLeft)/map->fSize > map->nColumns) {
                    nx = map->nColumns - 1;
                } else {
                    nx = (particle->fX - net->fLeft)/map->fSize;
                }
                
                if ((particle->fY - net->fBottom)/map->fSize < 0) {
                    ny = 0;
                } else if ((particle->fY - net->fBottom)/map->fSize > map->nRows) {
                    ny = map->nRows - 1;
                } else {
                    ny = (particle->fY - net->fBottom)/map->fSize;
                }
                
                map->pFlowMap[ny*map->nColumns+nx] += 1;
            }
            
            free(particle);
            particle  = NULL;
            
            // export flow maps
            if((fpM = fopen(pMPath, "a+"))==NULL) {
                printf("flow map file open error\n");
                fclose(fpM);
                return -1;
            } else {
                fprintf(fpM, "%d\t%d\n", map->nRows, map->nColumns);
                
                for(i=map->nRows-1; i>=0; i--) {
                    for(j=0;j<map->nColumns;j++) {
                        fprintf(fpM, "%d\t", map->pFlowMap[i*map->nColumns+j]);
                    }
                    fprintf(fpM, "\n");
                }
                fclose(fpM);
            }
            
            // export outlet discharges
            if((fpM = fopen(pOPath, "a+"))==NULL) {
                printf("outlet file open error\n");
                fclose(fpM);
                return -1;
            } else {
                fprintf(fpM, "%d\t%d\n", nT, nMP);
                fclose(fpM);
            }
        }
        
        finish_export = clock();
        fDuration_export += (double)(finish_export - start_export);
        
        nT = nT + nIT;
    }
    
    finish = clock();
    fDuration += (double)(finish - start);
    
    if (comm_sz == 1) {
        fDuration = fDuration - fDuration_export;
    }
    
    // export computational time cost
    printf("%d %d %f\n", comm_sz, runoff->nNumParticles, fDuration);
    if((fpM = fopen(pTPath, "a+"))==NULL) {
        printf("time cost file open error\n");
        fclose(fpM);
        return -1;
    } else {
        fprintf(fpM, "%d %d %f\n", comm_sz, runoff->nNumParticles, fDuration);
        fclose(fpM);
    }
    
    fpM = NULL;
    free(pOPath);
    pOPath = NULL;
    free(pMPath);
    pMPath = NULL;
    free(pTPath);
    pTPath = NULL;
    
    free(runoff->pRP_List);
    runoff->pRP_List = NULL;
    free(runoff);
    runoff = NULL;
    
    free(map->pFlowMap);
    map->pFlowMap = NULL;
    free(map);
    map = NULL;
    
    for (i=0;i<pPPMap->nNumPoints; i++) {
        free(pPPMap->pPPR_List[i].pfPPR);
        pPPMap->pPPR_List[i].pfPPR = NULL;
    }
    
    free(pPPMap->pPPR_List);
    pPPMap->pPPR_List = NULL;
    free(pPPMap);
    pPPMap = NULL;
    
    free(net->pLid_List);
    net->pLid_List = NULL;
    
    free(net->pBNid_List);
    net->pBNid_List = NULL;
    
    free(net->pSNid_List);
    net->pSNid_List = NULL;
    
    free(net->pNid_List);
    net->pNid_List = NULL;
    
    free(net->pLine_List);
    net->pLine_List = NULL;
    
    free(net->pNode_List);
    net->pNode_List = NULL;
    
    free(net);
    net = NULL;
    
    MPI_Finalize();
    
    return 0;
}
