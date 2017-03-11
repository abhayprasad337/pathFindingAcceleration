/*
ABSTRACT: This code is used to find the shortest path between two points
          in a 2-D plane.
*/


#include <iostream>
#include <string>
using namespace std;


int findNeighbor(int x, int y, int distToXY, int *neighborX, int*neighborY, int *dist, int length, int width, int* closedListX, int* closedListY, int NclosedList);
int findUnique(int * listX,int * listY, int *distXY, int openListcnt);
int removeElement(int *arrayList, int pos, int numElements);

int LENGTH = 6;
int WIDTH = 6;
int main()
{
    int Ncount;
    int neighborsX[8],neighborsY[8];
    int dists[8]={0,0,0,0,0,0,0,0};

    int closedListsX[10];///={4,5};
    int closedListsY[10];///={4,4};
    int closedListDists[10];

    int openListsX[100];
    int openListsY[100];
    int openListDists[100];

    int *closedListX = &closedListsX[0];
    int *closedListY = &closedListsY[0];
    int *closedListDist = &closedListDists[0];

    int *neighborX = &neighborsX[0];
    int *neighborY = &neighborsY[0];
    int *dist = &dists[0];

    int startX = 3;
    int startY = 3;

    int endX = 6;
    int endY = 6;

    int pathFound = 1;
    int openListCnt = 0;
    int closedListCnt = 0;

    Ncount = findNeighbor(startX, startY, 0, neighborX, neighborY, dist , LENGTH, WIDTH, closedListX, closedListY, closedListCnt);

    /// Add start node to closed list
    closedListsX[0] = startX;
    closedListsY[0] = startY;
    closedListCnt++;

    for (int i = 0;i<Ncount;i++){
        cout << neighborX[i] << "," << neighborY[i] << "-->" << dist[i] << endl;
    }


    /// Add first neighbors to closed list
    for (int addToClosedList = 0; addToClosedList < Ncount; addToClosedList++){
         closedListsX[closedListCnt] = neighborsX[addToClosedList];
         closedListsY[closedListCnt] = neighborsY[addToClosedList];
         closedListDists[closedListCnt++] = dists[addToClosedList];
    }


    /// Clear neighbors before finding new list of neighbors
    /*  This has been taken care as Ncount is returned from the
        findNeighbors function, care must be taken at the programmer's
        end to evaluate the neighbor array only till the value of Ncount */

    /// Add neighbors to openList
    for (int evaluateOpenList = 1; evaluateOpenList < closedListCnt ;evaluateOpenList++){
        Ncount = findNeighbor(closedListsX[evaluateOpenList],closedListsY[evaluateOpenList],closedListDists[evaluateOpenList],neighborX, neighborY,dist,LENGTH,WIDTH,closedListX,closedListY,closedListCnt);
        for (int addToOpenList = 0; addToOpenList < Ncount; addToOpenList++){
            openListsX[openListCnt] = neighborsX[addToOpenList];
            openListsY[openListCnt] = neighborsY[addToOpenList];
            openListDists[openListCnt++] = dist[addToOpenList];

        }
    }
//    for (int i = 0;i<openListCnt;i++)
//        cout << i<< ":" << openListsX[i] << "," << openListsY[i] << endl;
//
//     cout << "Eliminating duplicates now..." << endl;

    /// Find only unique elements of openList
    openListCnt = findUnique(&openListsX[0],&openListsY[0],&openListDists[0],openListCnt);

//    /// Test with dist below ***
//    cout << "Updated openList now is.." << endl;
//    for (int i = 0;i<openListCnt;i++)
//        cout << i<< ":" << openListsX[i] << "," << openListsY[i] << "->" << openListDists[i] << endl;


    return 0;
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
            //cout << i << " " << tmpOpenListCountX << endl;
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
