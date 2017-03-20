/*
ABSTRACT: This code is used to find the shortest path between two points
          in a 2-D plane.
*/


#include <iostream>
#include <string>
#include <stdlib.h>
#include <cstddef>
#include <fstream>
#include <limits>
#include <map>
#include <omp.h>
#include <sys/time.h>



using namespace std;


int findNeighbor(int x, int y, int distToXY, int *neighborX, int*neighborY, int *dist, int length, int width, int* closedListX, int* closedListY, int NclosedList);
int findUnique(int * listX,int * listY, int *distXY, int openListcnt);
int removeElement(int *arrayList, int pos, int numElements);
int findNeighborWithoutElimiation(int x, int y, int distToXY, int *neighborX, int*neighborY, int *dist, int length, int width);
int findShortestDijkstra(int* g, int N, int *shortestPath, int *pathDistance);

int LENGTH = 5;
int WIDTH = 5;

int startX = 3;
int startY = 3;

int endX = 1;
int endY = 2;

typedef struct genPath{
    int pXval;
    int pYval;
    int x_neighbors[8];
    int y_neighbors[8];
    int dists[8];
    int numNeighbors;
    struct genPath * next;
}genPath;

typedef struct backTrack{
    int x;
    int y;
    int index;
    backTrack * next;
}backTrack_t;


int listLength (genPath* head);
int listLengthbt (backTrack_t* head);
int insertDataIntoLinkedList (genPath** head, int xVal, int yVal, int *neighborX, int * neighborY, int *dist, int neighborCount);
void reverseLinkedList (genPath** head);
void reverseLinkedListbt (backTrack_t** head);
void printGraph(int * graph,int xlim, int ylim);
int findIndex (backTrack_t* headBT,int x, int y);
void indexXY(genPath * head,backTrack_t **headBacktrack);
int findMinLoc(int dist[],bool pri[],int N);
int findXYfromIndex (backTrack_t* headBT, int index, int* Xval, int* Yval);
void printPath (int startLocX, int startLocY, int endLocX, int endLocY, int *pathX, int *pathY, int pathLength, int terrainLength, int terrainWidth);


int main()
{
    int Ncount;
    int neighborsX[8],neighborsY[8];
    int tmpNeighborsX[8], tmpNeighborsY[8];
    int dists[8]={0,0,0,0,0,0,0,0};
    int tmpDists[8] = {0,0,0,0,0,0,0,0};

    int closedListsX[1000];///={4,5};
    int closedListsY[1000];///={4,4};
    int closedListDists[1000];

    int openListsX[1000];
    int openListsY[1000];
    int openListDists[1000];

    int tmpListX[1000];
    int tmpListY[1000];
    int tmpListDists[1000];

    int *closedListX = &closedListsX[0];
    int *closedListY = &closedListsY[0];
    //int *closedListDist = &closedListDists[0];

    int *neighborX = &neighborsX[0];
    int *neighborY = &neighborsY[0];
    int *dist = &dists[0];

    int *tmpNeighborX = &tmpNeighborsX[0];
    int *tmpNeighborY = &tmpNeighborsY[0];
    int *tmpDist = &tmpDists[0];

    int pathFound = 0;
    int openListCnt = 0;
    int closedListCnt = 0;
    int mopenListCnt=0;
    int terminate = 0; /// This is used to indicate that the destination was found.
    int iterations =0; /// This is used to monitor the number of iterations required to reach destination.
    int tmpListCnt;

    /** Initialize genPath linkedList **/
    genPath *head = (genPath*) malloc(sizeof(genPath));
    if (head == NULL){
        return 1;
    }

    /* Add start node to closed list */
    closedListsX[0] = startX;
    closedListsY[0] = startY;
    closedListCnt++;

    /** Find neighbors to start pathfinding **/
    /**
     * This is done to load the initialize the loop for
     * building the path. This is the nature of the algorithm
     * that has been designed. A better algorithm for this
     * should be implemented which does not handle such cases
     * this bad.
    */

    Ncount = findNeighbor(startX, startY, 0, neighborX, neighborY, dist , LENGTH, WIDTH, closedListX, closedListY, closedListCnt);

    /// Manually insert into linkedlist for startnode
    head->next = NULL;
    head->pXval = startX;
    head->pYval = startY;
    head->numNeighbors = Ncount;
    memcpy(head->dists,dist,sizeof(int)*Ncount);
    memcpy(head->x_neighbors,neighborX,sizeof(int)*Ncount);
    memcpy(head->y_neighbors,neighborY,sizeof(int)*Ncount);

    /** Add the present neighbors of start node to closed list **/
    /** This is due to the nature of the algorithm. A much cleaner
     *  code could be written if the algorithm can be modified.
    */
    for (int addToClosedList = 0; addToClosedList < Ncount; addToClosedList++){
        closedListsX[closedListCnt] = neighborsX[addToClosedList];
        closedListsY[closedListCnt] = neighborsY[addToClosedList];
        closedListDists[closedListCnt++] = dists[addToClosedList];
    }

    /** Add the present neighbors of start node to open list **/
    /** This is due to the nature of the algorithm. A much cleaner
     *  code could be written if the algorithm can be modified.
    */
    for (int addToOpenList = 0; addToOpenList < Ncount; addToOpenList++){
        openListsX[openListCnt] = neighborsX[addToOpenList];
        openListsY[openListCnt] = neighborsY[addToOpenList];
        openListDists[openListCnt++] = dist[addToOpenList];
    }

    /**
      * Find unique elements on the openList.
      * Is this useful at all??? Clean up if necessary
      * - According to the review, this is done to sort of
      *   initialize mopenListcnt/
      */
    mopenListCnt = findUnique(&openListsX[0],&openListsY[0],&openListDists[0],openListCnt);


    /**
     * This loop is used to build the path till the destination is reached. A linkedlist node
     * for each element is created till the destination is reached.
    */

    while (!pathFound){

        cout << "Iteration = " << iterations++ << endl;

        int neighborCount = 0;/// Keeps track of number of valid neighbors found.
        int tmpNeighborCount = 0;/// Load tmpNeighbors for linkedList backtrack.

        tmpListCnt = 0;
        for (int i=0; i<mopenListCnt; i++){

            cout  << openListsX[i] << "," << openListsY[i] << endl;
            neighborCount = findNeighbor(openListsX[i], openListsY[i], 0, neighborX, neighborY, dist , LENGTH, WIDTH, closedListX, closedListY, closedListCnt);

            // Declare another variable and put into linkedList.
            tmpNeighborCount = findNeighborWithoutElimiation(openListsX[i], openListsY[i], 0, tmpNeighborX, tmpNeighborY, tmpDist, LENGTH, WIDTH);


            if(!insertDataIntoLinkedList (&head,openListsX[i],openListsY[i],&tmpNeighborsX[0],&tmpNeighborsY[0],&tmpDist[0],tmpNeighborCount)){
                cout << "ERROR: Failed to add item to LinkedList" << endl;
            }

            /**
             * Temporarily stores the list of all neighbors.
             * This can contain duplicate values.
            */
            for (int addToTmpList = 0; addToTmpList < neighborCount; addToTmpList++){
                tmpListX[tmpListCnt] = neighborsX[addToTmpList];
                tmpListY[tmpListCnt] = neighborsY[addToTmpList];
                tmpListDists[tmpListCnt++] = dist[addToTmpList];
            }



            /**
             * To check if the destination node is reached.
             * This is done in the (n+1)th iteration as there
             * are chances of us missing out on shorter paths
             * in the n-th level.
            */
            if ((openListsX[i]==endX) && (openListsY[i]==endY)){
                cout << "caught you!" << endl;
                terminate = 1;
                break;
            }
        }



        if (!terminate){

            /**
             * Find unique elements in the tmpList. tmpList
             * contains neighbors of every pixel in the nth
             * level of iteration. In order to find the next
             * few nodes to evaluate, we need to find the unique
             * of these.
            */
            mopenListCnt = findUnique(&tmpListX[0],&tmpListY[0],&tmpListDists[0],tmpListCnt);

            /**
              * Add all unique elements to the openList.
              * Note that the openList index starts from 0 every time.
              * So this completely refreshed for each iteration.
            */
            for (int addToOpenList = 0; addToOpenList < mopenListCnt; addToOpenList++){
                openListsX[addToOpenList] = tmpListX[addToOpenList];
                openListsY[addToOpenList] = tmpListY[addToOpenList];
                openListDists[addToOpenList] = tmpListDists[addToOpenList];
                //cout << "On iteration " << iterations << " Values are = " << openListsX[addToOpenList] << "," << openListsY[addToOpenList] << "-->" << openListDists[addToOpenList] << endl;
            }

            /**
             * Add elements to the closed list. This list would help in
             * eliminating the points that have already been evaluated.
             * Notice that this is being done in the (n+1)th iteration.
            */
            for (int addToClosedList = 0; addToClosedList < mopenListCnt; addToClosedList++){
                closedListsX[closedListCnt] = tmpListX[addToClosedList];
                closedListsY[closedListCnt] = tmpListY[addToClosedList];
                closedListDists[closedListCnt++] = tmpListDists[addToClosedList];
            }


        } else {
            pathFound = 1;
        }

       cout << "List length = " << listLength (head)  << endl;
    }



    /** Reverse the order of linkedlist */
    reverseLinkedList (&head);
    int totalNodes = listLength(head);


    cout << "Traced end point. " << "Found " <<  totalNodes << " nodes. Starting to backtrack.." << endl;

    /**
     * At this point of the code, all nodes from the start to the destination have been
     * traversed. Backtracking must be done to build the distance matrix.
    */

    const int xgraphWidth = totalNodes;
    const int ygraphWidth = totalNodes;
    int graph[xgraphWidth][ygraphWidth]={0};

    backTrack_t * headBacktrack = (backTrack_t*) malloc(sizeof(backTrack_t));
    headBacktrack->next = NULL;
    headBacktrack->x = head->pXval;
    headBacktrack->y = head->pYval;
    headBacktrack->index = 0;

    /// Index all X and Y locations
    indexXY(head,&headBacktrack);
    reverseLinkedListbt(&headBacktrack);
    cout << "Length of indexed LL " << listLengthbt(headBacktrack) << endl;

    genPath * current = head;
    int index;
    int counter = 0;

    while (current!=NULL){

        //cout << current->pXval << "," << current->pYval << endl;

        cout << "Neighbors of " << current->pXval << "," << current->pYval << " is.." << endl;

        for (int thisIndex = 0; thisIndex<current->numNeighbors; thisIndex++){
                cout << current->x_neighbors[thisIndex] << "," << current->y_neighbors[thisIndex] << endl;
                index = findIndex(headBacktrack,current->x_neighbors[thisIndex],current->y_neighbors[thisIndex]);
                if (index >= 0){
                    graph[index][counter] = current->dists[thisIndex];
                } else {
                    graph[index] [counter]  = 0;
                    //cout << "Could not find the index." << endl;
                }
        }

        cout << "counter now is = " << counter << endl;
        counter ++;
        current = current->next;
    }

    //printGraph(&graph[0][0],totalNodes, totalNodes);
    for (int i=0;i<totalNodes;i++){
        for (int j=0;j<totalNodes;j++){
            cout << graph[i][j] << " ";
        }
        cout << endl;
    }

    int pathDistance = 0;
    int shortestPath[totalNodes]={-1};
    int Nvals;
    Nvals = findShortestDijkstra(&graph[0][0], totalNodes, &shortestPath[0], &pathDistance);

//    cout << "Path is = " << endl;
//    for (int i = 0; i < Nvals; i++){
//        cout << shortestPath [i] << endl;
//    }


    int Xval[Nvals],Yval[Nvals];
    for (int i = 0; i < Nvals; i++){
        Xval[i] = -1;
        Yval[i] = -1;
        if(findXYfromIndex (headBacktrack, shortestPath[i],&Xval[i],&Yval[i])){
                cout << shortestPath[i] << " ---> " << Xval[i] << "," << Yval[i] << endl;
        } else {
                cout << "Point not found in linked list" << endl;
        }

    }
    cout << "Distance to reach goal = " << pathDistance << endl;

    cout << "Printing result to file..." << endl;
    printPath (startX, startY, endX, endY, &Xval[0], &Yval[0], Nvals, LENGTH, WIDTH);

    return 0;
}
void printPath (int startLocX, int startLocY, int endLocX, int endLocY, int *pathX, int *pathY, int pathLength, int terrainLength, int terrainWidth){
       int terrain[terrainLength][terrainWidth]={0};
       int Xval,Yval;
       ///terrain[terrainLength-1][terrainWidth-1]=0;

       for (int i = 0;i < pathLength; i++){
            Xval = pathX[i];
            Yval = pathY[i];
            terrain[Xval][Yval] = 127;
       }

       terrain[startLocX][startLocY] = 255;
       terrain[endLocX][endLocY] = 255;

       for (int i=0;i<terrainLength;i++){
            for (int j=0; j<terrainWidth;j++){
                cout << terrain[i][j] << ",";
            }
            cout << endl;
       }
       /// Print into file -> path.txt
        ofstream myfile;
        myfile.open ("path.txt");

        for (int i=0;i<terrainLength;i++){
            for (int j=0;j<terrainWidth;j++){
                myfile << terrain[i][j] << " ";
            }
            myfile << endl;
        }

        myfile.close();
        return;
}

int findXYfromIndex (backTrack_t* headBT, int index, int* Xval, int* Yval){
    while(headBT!=NULL){
        if(headBT->index==index){
               *Xval = headBT->x;
               *Yval = headBT->y;
               return 1;
        }
        headBT = headBT->next;
    }
    return -1;
}

void indexXY(genPath * head,backTrack_t **headBacktrack){
    int incrementer = 1;
    head = head->next;
    while(head->next!=NULL){
            backTrack_t* newNode;
            newNode = (backTrack_t*) malloc(sizeof(backTrack_t));
            newNode->next = *headBacktrack;
            newNode->x = head->pXval;
            newNode->y = head->pYval;
            newNode->index = incrementer++;
            *headBacktrack = newNode;
            head = head->next;
    }

    backTrack_t* newNode;
    newNode = (backTrack_t*) malloc(sizeof(backTrack_t));
    newNode->next = *headBacktrack;
    newNode->x = head->pXval;
    newNode->y = head->pYval;
    newNode->index = incrementer++;
    *headBacktrack = newNode;
}

int findIndex (backTrack_t* headBT,int x, int y){
    while(headBT!=NULL){
        if(headBT->x==x){
                if (headBT->y==y){
                        return headBT->index;
                }
        }
        headBT = headBT->next;
    }
    return -1;
}

int listLength (genPath* head){
    genPath * current = head;
    int count =0;
    while (current != NULL){
        count ++;
        //cout << "Monitor HEAD: " << current->pXval << "," << current->pYval << endl;
        current = current -> next;
    }
    return count;
}

int listLengthbt (backTrack_t* head){
    backTrack_t * current = head;
    int count =0;
    while (current != NULL){
        count ++;
        cout << "Monitor: " << current->x << "," << current->y << " -> "<< current->index << endl;
        current = current -> next;
    }
    return count;
}

int insertDataIntoLinkedList (genPath** head, int xVal, int yVal, int *neighborX, int * neighborY, int *dist, int neighborCount){
    genPath* newNode;
    newNode = (genPath*) malloc(sizeof(genPath));

    if (!newNode){
        cout << "Memory error!" << endl;
        return 0;
    }

    newNode->next = *head;
    newNode->pXval = xVal;
    newNode->pYval = yVal;
    newNode->numNeighbors = neighborCount;
    memcpy(newNode->dists,dist,sizeof(int)*neighborCount);
    memcpy(newNode->x_neighbors,neighborX,sizeof(int)*neighborCount);
    memcpy(newNode->y_neighbors,neighborY,sizeof(int)*neighborCount);
    *head = newNode;
    return 1;
}

void reverseLinkedList (genPath** head){
    genPath* currentNode = *head;
    genPath* nextNode = NULL;
    genPath* prevNode = NULL;

    while(currentNode!=NULL){
		nextNode = currentNode->next;
		currentNode->next = prevNode;
		prevNode = currentNode;
		currentNode = nextNode;
	}
	*head = prevNode;
    return;
}

void reverseLinkedListbt (backTrack_t** head){
    backTrack_t* currentNode = *head;
    backTrack_t* nextNode = NULL;
    backTrack_t* prevNode = NULL;

    while(currentNode!=NULL){
		nextNode = currentNode->next;
		currentNode->next = prevNode;
		prevNode = currentNode;
		currentNode = nextNode;
	}
	*head = prevNode;
    return;
}

int findUnique(int * listX,int * listY, int *distXY, int openListcnt){
    int x_check,y_check;
    int duplicate=0;
    int checkedOff[openListcnt]={0};
    for (int i=0;i<openListcnt;i++){
        if (checkedOff[i]==0){
            x_check = listX[i];
            y_check = listY[i];
            for (int j=i+1;j<openListcnt;j++){
                if (x_check == listX[j] && y_check == listY[j] && checkedOff[j]==0){
                    checkedOff[j] = 1;
                    duplicate++;
                }
            }
        }
    }
    cout << "number of duplicates found = " << duplicate << endl;

    int tmpOpenListCountX=openListcnt;
    int tmpOpenListCountY=openListcnt;
    int tmpDistCount = openListcnt;
    int delCount = 0;
    for (int i=0;i<openListcnt;i++){
       if (checkedOff[i]==1){
           // cout << "DUBLICATE " << i << " " << listX[tmpOpenListCountX] << "," << listY[tmpOpenListCountX] << endl;
            tmpOpenListCountX = removeElement(listX,i-delCount,tmpOpenListCountX);
            tmpOpenListCountY = removeElement(listY,i-delCount,tmpOpenListCountY);
            tmpDistCount = removeElement(distXY,i-delCount,tmpDistCount);
            delCount++;
        }
    }

    /// Full proof the function
    if (tmpOpenListCountX!=tmpOpenListCountY){
        tmpOpenListCountX = -1;
    }
    /// This is the number of elements after eliminating the
    /// duplicate ones.
    return tmpOpenListCountX;
}

 /* ABSTRACT: This function is used to remove an array element.
  *           It is similar to the pop() function.
  *
  * INPUTS:
  *     arrayList = Pointer to the first element of the array.
  *     pos = Position of element to remove. NOTE: This is
  *           a 0 indexed value. To remove the 3rd element
  *           from the array, pos should be 2.
  *     numElements = Total number of elements of the array.
  * OUTPUTS:
  *    rnumElemnts = Number of elements after pop().
  */
int removeElement(int *arrayList, int pos, int numElements){
   memmove(&arrayList[pos],&arrayList[pos+1],sizeof(int)*(numElements-pos-1));
   int rnumElements = numElements -1;
   return rnumElements;
}

/* ABSTRACT: This function is used to find the neighbors of
 *           a given pixel in the grid. Neighbor pixels found
             in closed list are omitted.

 * INPUTS: x - Pixel x location in the grid
           y - Pixel y location in the grid
           neighborX - An empty 8 element array to store the X indices neighboring pixels.
           neighborY - An empty 8 element array to store the Y indices of neighboring pixels.
           dist - An 8 element empty array.
           length - Length of the entire grid. This is used to eliminate neighbor pixels beyond the grid size.
           width - Width of the entire grid. This is used to eliminate neighbor pixels beyond the grid size.
           closedListX - X indices of the elements of the closed list.
           closedListY - Y indices of the elements of the closed list.
           NclosedList - Number of element sin the closed list.

  OUTPUTS: Ncount - Number of valid neighbors for the given pixel.
*/
int findNeighbor(int x, int y, int distToXY, int *neighborX, int*neighborY, int *dist, int length, int width, int* closedListX, int* closedListY, int NclosedList){

    int validX[3],validY[3];
    int xVals[3] = {x-1,x,x+1};
    int yVals[3] = {y-1,y,y+1};


    if (xVals[0] < 0 || xVals[0] > width){
        validX[0] = 0;
    } else {
        validX[0] = 1;
    }

    if (xVals[2] < 0 || xVals[2] > width){
        validX[2] = 0;
    } else {
        validX[2] = 1;
    }

    if (yVals[0] < 0 || yVals[0] > length){
        validY[0] = 0;
    } else {
        validY[0] = 1;
    }

    if (yVals[2] < 0 || yVals[2] > length){
        validY[2] = 0;
    } else {
        validY[2] = 1;
    }

    validX[1] = 1;
    validY[1] = 1;

    int Ncount=0;
    int dontConsider=0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            if (validX[i] && validY[j] && !(i==1&&j==1)){
            dontConsider = 0;
              for (int k=0 ;k < NclosedList; k++){
                if (xVals[i] == closedListX[k] && yVals[j] == closedListY[k]){
                    dontConsider = 1;
                }
              }
              if (!dontConsider){
                neighborX[Ncount] = xVals[i];
                neighborY[Ncount] = yVals[j];
                if ((i==0 && j==0) || (i==0&&j==2) || (i==2&&j==0) || (i==2&&j==2)){
                    dist[Ncount] = distToXY + 14;
                } else {
                    dist[Ncount] = distToXY + 10;
                }
                Ncount++;
              }
            }
        }
    }

    return Ncount;
}

int findNeighborWithoutElimiation(int x, int y, int distToXY, int *neighborX, int*neighborY, int *dist, int length, int width){

    int validX[3],validY[3];
    int xVals[3] = {x-1,x,x+1};
    int yVals[3] = {y-1,y,y+1};


    if (xVals[0] < 0 || xVals[0] > width){
        validX[0] = 0;
    } else {
        validX[0] = 1;
    }

    if (xVals[2] < 0 || xVals[2] > width){
        validX[2] = 0;
    } else {
        validX[2] = 1;
    }

    if (yVals[0] < 0 || yVals[0] > length){
        validY[0] = 0;
    } else {
        validY[0] = 1;
    }

    if (yVals[2] < 0 || yVals[2] > length){
        validY[2] = 0;
    } else {
        validY[2] = 1;
    }

    validX[1] = 1;
    validY[1] = 1;

    int Ncount=0;
//    int dontConsider=0;
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
            if (validX[i] && validY[j] && !(i==1&&j==1)){
            //dontConsider = 0;
//              for (int k=0 ;k < NclosedList; k++){
//                if (xVals[i] == closedListX[k] && yVals[j] == closedListY[k]){
//                    dontConsider = 1;
//                }
//              }
              //if (!dontConsider){
                neighborX[Ncount] = xVals[i];
                neighborY[Ncount] = yVals[j];
                if ((i==0 && j==0) || (i==0&&j==2) || (i==2&&j==0) || (i==2&&j==2)){
                    dist[Ncount] = distToXY + 14;
                } else {
                    dist[Ncount] = distToXY + 10;
                }
                Ncount++;
              //}
            }
        }
    }

    return Ncount;
}

int findMinLoc(int dist[],bool pri[],int N){
    int minval = INT_MAX;
    int minloc = 0;;
    for (int i = 0; i< N; i++){
        if(!pri[i] && dist[i]<minval){
            minval = dist[i];
            minloc = i;
        }
    }
    return minloc;
}

int findShortestDijkstra(int* tmp_g, int N, int *shortestPath, int *pathDistance){

    int i,j;
    int dist[N];
    bool pri[N];
    int path[N];

    int g[N][N]={0};
    struct timeval start_time, stop_time, elapsed_time;  // timers


    for (int i=0;i<N;i++){
        for (int j=0;j<N;j++){
            g[i][j] = *((tmp_g+i*N) + j);
        }
    }
//
//    for (int i=0;i<N;i++){
//        for (int j=0;j<N;j++){
//            cout << g[i][j] << ",";
//        }
//        cout << endl;
//    }

    /// Initialize all values in the queue infinity.
    for (i = 0 ; i < N ; i++){
        dist[i] = INT_MAX; /// Infinity
        pri[i] = false;
        path[i] = 0;
    }

    /// Make distance from source to source as 0
    dist[0] = 0;
    path[0] = -1; /// Make source to be -1
    int u;


   /// Path finding..
   gettimeofday(&start_time,NULL); // Unix timer
    for (i = 0 ; i < N-1 ; i++){
        u = findMinLoc(dist,pri,N);
        pri[u] = true;

        for (j = 0; j < N ;j++){
                if (g[u][j] && g[u][j] + dist[u] < dist[j] && dist[u]!=INT_MAX && !pri[j]){
                    /// update queue
                    dist[j] = g[u][j] + dist[u];
                    path[j] = u;
                }
        }
    }

    gettimeofday(&stop_time,NULL);
//    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract
    printf("Total time was %f seconds.\n", elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0);

    for (i = 0 ; i < N ; i++){
        cout << dist[i] <<endl;
    }
    *pathDistance = dist[N-1];

//    cout << "stored in path format --> " << endl;
//    for (i = 0 ; i < N ; i++){
//        cout << path[i] <<endl;
//    }


    /// Backtrack to find the right path.
    int traceBack = N-1;
    int finalPath[N];

    /// Initialize array to -1.
    for (i = 0 ; i < N ; i++){
        finalPath[i] = -1;
    }

    finalPath[0] = traceBack;
    i=1;

    /// Track path until source is found.
    while(traceBack > 0 ){
        finalPath[i] = path[traceBack];
        traceBack =  finalPath[i];
        i++;
    }

    int Nvals=0;
    /// Print the path.
    cout << "path successfully found --> " << endl;
    for (i = 0; i < N; i++){
        if (finalPath[i]!=-1){
            *shortestPath++ = finalPath[i];
            Nvals = i;
        }
        cout << finalPath[i] <<endl;
    }
    return ++Nvals;
}
