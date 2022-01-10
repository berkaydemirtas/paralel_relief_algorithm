//Student Name: Berkay Demirtaş
//Student Number: 2017400234
//Compile Status: Compiling
//Program Status: Working

#include <iostream>
#include <mpi.h>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include "mpi.h"
#include <fstream>
#include <vector>
#include <bits/stdc++.h>
#include <algorithm>

using namespace std;

int N; // number of instance
int A; // number of features
int M; //iteration count
int T; // top features number
int P;  // total process number


//it is a simple function calculates manhattan distance between 2 vectors
double calculateManhattanDistance(vector<double> v1, vector<double> v2){
    double sum=0;
    for(int i=0; i< A; i++){
        sum+=abs(v1[i]-v2[i]);
    }
    return sum;
}


//takes pref and returnValues vectors and selected int as an argument. Selected corresponds to id of the pivot instance
//return values used to return ids of nearest hit and nearest miss
//for each instance in pref vector, it calculates manhattan distance with pivot instance by using manhattanDistance function
//if class is same and value is smaller than smallest value until this time, then min value is changed
//and same prosedure for different class
//at the end finded ids are returned with an vector.
vector<int> calculateNearestHitandMiss(vector<double> pref,vector<int> returnValues, int selected){
    double sameClassMinValue=INFINITY;
    double sameClassMinInstance=-1;
    double diffClassMinValue=INFINITY;
    double diffClassMinInstance=-1;
    vector<double> selectedFeatures;
    int selectedClass=pref[(selected+1)*(A+1)-1];
    for(int i=selected*(A+1);i<(selected+1)*(A+1)-1;i++){
        selectedFeatures.push_back(pref[i]);
    }
    for(int i=0;i<N/(P-1);i++){

        vector<double> selectedFeatures2;
        if(i!=selected){

            int nextClass=pref[(i+1)*(A+1)-1];
            for(int j=i*(A+1);j<(i+1)*(A+1)-1;j++){
                selectedFeatures2.push_back(pref[j]);
            }

            double result = calculateManhattanDistance(selectedFeatures,selectedFeatures2);

            if(selectedClass==nextClass){
                if(result<sameClassMinValue){
                    sameClassMinValue=result;
                    sameClassMinInstance=i;
                }
            }
            else{
                if(result<diffClassMinValue){
                    diffClassMinValue=result;
                    diffClassMinInstance=i;
                }
            }
        }
    }
    returnValues.push_back(sameClassMinInstance);
    returnValues.push_back(diffClassMinInstance);

    return returnValues;
}

// takes pref vector and feature num int as arguments and find maximum and minimum of
// values in pref vector with featureNum values.
//returns difference between them.
double MaxMin(vector<double> pref, int featureNum){
    double max=-INFINITY;
    double min = INFINITY;
    for(int i=0;i<N/(P-1);i++){
        double currentFeature = pref[i*(A+1)+featureNum];
        if(currentFeature<min)
            min = currentFeature;
        if(currentFeature>max)
            max=currentFeature;
    }
    return max-min;
}


// takes W vector which contains relief weights all with 0 and updates it.
//first for is iterates iteration count times and second for calculates each W[a] value for
//each iteration.
//this function takes nearestHit and nearestMiss ids from calculateNearestHitandMiss function and
// also uses the MaxMin function.
void relief(vector<double> &W,vector<double> pref){

    for(int i=0;i<M;i++){

        vector<int> returnValues;
        returnValues=calculateNearestHitandMiss(pref,returnValues,i%(N*(A+1)/(P-1)));
        int nearestHit=returnValues[0];
        int nearestMiss=returnValues[1];
        for(int a=0;a<A;a++) {

            double currentFeatureR = pref[i%(N*(A+1)/(P-1)) * (A + 1) + a];
            double currentFeatureH = pref[nearestHit * (A + 1) + a];
            double currentFeatureM = pref[nearestMiss * (A + 1) + a];
            W[a]= W[a]- (abs(currentFeatureR-currentFeatureH)/MaxMin(pref,a))/M + (abs(currentFeatureR-currentFeatureM)/MaxMin(pref,a))/M ;
        }


    }
}

int main(int argc, char* argv[])
{
    int rank; // rank of the current processor

    MPI_Init(&argc, &argv); //mpi starts
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // gets the rank of the current processor
    MPI_Comm_size(MPI_COMM_WORLD, &P); // gets the total number of processors


    // argv[1] is the input file
    // and this identifiers are explained above
    std::fstream fileStream;
    fileStream.open(argv[1]);
    fileStream >> P;
    fileStream >> N;
    fileStream >> A;
    fileStream >> M;
    fileStream >> T;

    double arr[(P*N*(A+1))/(P-1)]; //array used in scatter
    double pref[N*(A+1)/(P-1)]; // local disk on processors
    int masterIndexes[T*P];// disk of master for mpi_gather
    int indexes2[T]; // used just for ease of implementing

    // If it's master processor, reads from input file
    if(rank==0){
        for(int i=0;i<N*(A+1)/(P-1);i++)
            arr[i]=0;
        for(int i=N*(A+1)/(P-1);i<(P*N*(A+1))/(P-1);i++){
            fileStream >> arr[i];
        }

        // sends M and T values to all slaves
        for(int i=1;i<P;i++){

            MPI_Send(
                    /* data         = */ &M,
                    /* count        = */ 1,
                    /* datatype     = */ MPI_INT,
                    /* destination  = */ i,
                    /* tag          = */ 0,
                    /* communicator = */ MPI_COMM_WORLD);

            MPI_Send(
                    /* data         = */ &T,
                    /* count        = */ 1,
                    /* datatype     = */ MPI_INT,
                    /* destination  = */ i,
                    /* tag          = */ 1,
                    /* communicator = */ MPI_COMM_WORLD);
        }
    }

    // sends data from root array arr to pref array on each processor
    MPI_Scatter(arr,N*(A+1)/(P-1),MPI_DOUBLE,pref,N*(A+1)/(P-1),MPI_DOUBLE,0,MPI_COMM_WORLD);

    int masterSignal = 1;
    while(masterSignal){

        //if it is a slave processor it recieves M and T values
        if(rank!= 0){
            MPI_Recv(
                    /* data         = */ &M,
                    /* count        = */ 1,
                    /* datatype     = */ MPI_INT,
                    /* source       = */ 0,
                    /* tag          = */ 0,
                    /* communicator = */ MPI_COMM_WORLD,
                    /* status       = */ MPI_STATUS_IGNORE);

            MPI_Recv(
                    /* data         = */ &T,
                    /* count        = */ 1,
                    /* datatype     = */ MPI_INT,
                    /* source       = */ 0,
                    /* tag          = */ 1,
                    /* communicator = */ MPI_COMM_WORLD,
                    /* status       = */ MPI_STATUS_IGNORE);

            M=M;
            T=T;

            //vector to use in mpi_gather
            vector<int> indexes;
            int i = 0;
            vector<double> pref2;
            for(; i<N*(A+1)/(P-1);i++) {
                pref2.push_back(pref[i]);
            }

            // W vector for relief algorithm
              vector<double> W;
              for(int j=0;j<A;j++){
                  W.push_back(0);
              }

              //calls relief function
              relief(W,pref2);

              //used for ease of coding
              vector<double> W2;

            cout<<"Slave P"+to_string(rank)+" : ";

              for(int j=0;j<W.size();j++){
                  W2.push_back(W[j]);
              }
              sort(W.begin(), W.end(), greater<double>());
              // sorts updated W vector and fills indexes vector
              for(int j=0;j<T;j++){
                  for(int k=0;k<A;k++) {
                      if (W2[k] == W[j])
                          indexes.push_back(k);
                  }
              }

                sort(indexes.begin(), indexes.end());
                //prints the features ids of current slave
                for(int j=0;j<indexes.size();j++) {
                    indexes2[j]=indexes[j];
                    if(j!=indexes.size()-1)
                    cout << indexes[j] << " ";
                }

                cout<<indexes[indexes.size()-1];

                cout<<endl;

        }

        //gather all feature ids from slaves to master
        MPI_Gather(indexes2, T, MPI_INT, masterIndexes, T, MPI_INT, 0,
                   MPI_COMM_WORLD);

        //in this part master removes dublicates from masterIndexes array and print indexes in ascending order
        if(rank==0){
            set<int> s1;
            vector<int> printer;
            cout<<"Master P0 : ";
            for(int i=T;i<P*T;i++){
                s1.insert(masterIndexes[i]);
            }
            set<int, greater<int> >::iterator itr;
            for (itr = s1.begin(); itr != s1.end(); itr++)
            {
                printer.push_back(*itr);
            }

            for(int i=0;i<printer.size()-1;i++)
                cout<<printer[i]<<" ";

            cout<<printer[printer.size()-1];

            cout<<endl;
            masterSignal=0;
        }

        MPI_Bcast(&masterSignal, 1, MPI_INT, 0, MPI_COMM_WORLD); // broadcast


    }


    // ****************************************** //


    MPI_Barrier(MPI_COMM_WORLD); // synchronizing processes
    MPI_Finalize();

    return 0;
}





