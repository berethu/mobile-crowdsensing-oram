#ifndef HASHFUNC_H
#define HASHFUNC_H


#include <string>
using namespace std;
//string arr[80000][3];
class hashFunc
{
    private:
        static const int partitionSize= 10000;

        struct item{
        string location;
        int speed;
        item* next;
        };

        item* HashTable[partitionSize];
    public:
        int pSize;
        hashFunc();
        int Hashfunc(string key);
        void AddItem(string location , int speed);
        int NumberOfItemsInIndex(int index);
        void PrintPartition();
        int SumAndPrintItemsInIndex(int index);
        int getData(string filelocation, int rowA, int colA, string arr[][3]);
        
        //string arr[80000][3]; //44,589,366



};


#endif // HASHFUNC_H
