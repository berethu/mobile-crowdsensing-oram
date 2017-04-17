//For aggregation. It adds all the data belonging to thesame partition together and pass it to SP
#include <cstdlib>
#include <iostream>
#include <string>
#include "hashFunc.h"
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <unistd.h>

using namespace std;

hashFunc::hashFunc()
{
    pSize = partitionSize;
    for(int i = 0; i<partitionSize; i++)
    {
        HashTable[i] = new item;
        HashTable[i]->location = "empty";
        HashTable[i]->speed = 0;
        HashTable[i]->next = NULL;
    }
}

int hashFunc::getData(string filelocation, int rowA, int colA, string arr[][3])
{
        string lineA;
     
    ifstream fileN;
    string x;
    //={{0}};
   // cout << "'Input file location" << endl;
    
     fileN.open(filelocation);

     //Error check
     if(fileN.fail())
     {
         cerr<<"file is bad";
         exit(1);
     }


     //Reading data
     cout<<"Reading data \n"<<endl;
     while (fileN.good())
     {
         while(getline(fileN, lineA))
         {
         istringstream streamA(lineA);
         colA = 0;
         while(streamA >>x) {
            arr[rowA][colA] = x;
            colA++;
         }
         rowA++;

         }
     }

     //Display data
     cout<<"# of rows ====> "<<rowA<<endl;
     cout<<"# of cols ===> "<<colA<<endl;
     cout<<" "<<endl;
     // for (int i =0; i<rowA; i++)
     // {
     //     for(int j =0; j<colA; j++)
     //     {
     //        cout<<left <<setw(6)<<arr[i][j]<<" ";
     //     }
     //     cout<<endl;
     // }
double percentage = 0;
cout<<"Reading data"<<endl;
   //  for (int i =0; i<rowA; i++)
     //{
   //     percentage =(i+1/rowA) * 10;// * 100);
        // for(int j =0; j<colA; j++)
        // {
        //    cout<<i+1 <<" of "<<rowA<<" percentage completion "<<percentage<<" %";
      //   }
         cout<<endl;
     //}
return rowA;
}



void hashFunc::AddItem(string location, int speed)
{
    int index = Hashfunc(location);
    if( HashTable[index]->location == "empty")
    {
        HashTable[index]->location = location;
        HashTable[index]->speed = speed;
    }
    else {
        item* Ptr= HashTable[index];
        item* n = new item;
        n->location = location; //+ HashTable[index]->location;
        n->speed = speed; //+ HashTable[index]->speed; //This does the aggregation
        n->next = NULL;
        while (Ptr->next !=NULL)
        {
            Ptr = Ptr->next;
        }
        Ptr->next = n;
    }
}

int hashFunc::NumberOfItemsInIndex(int index)
{
int count = 0;

if (HashTable[index]->location == "empty")
{
    return count;
}
else{
    count++;
    item* Ptr = HashTable[index];
    while(Ptr->next != NULL)
    {
        count++;
        Ptr = Ptr->next;
    }
}
    return count;
}

void hashFunc::PrintPartition()
{
    int number;
    for(int i=0;i<partitionSize; i++)
    {
        number = NumberOfItemsInIndex(i);
        cout<<"-------------\n";
        cout<<"index = "<<i<<endl;
        cout<<HashTable[i]->location<<endl;
        cout<<HashTable[i]->speed<<endl;
        cout<<"# of items = "<<number<<endl;
        cout<<"-----------------\n";
    }
}


int hashFunc::Hashfunc(string key)
{

    int hashnum =0;
    int index;

    for (int i = 0; i<key.length(); i++)
    {
        hashnum = hashnum + (int)key[i];
        cout<<"hash["<<i<<"] = "<< hashnum<<endl;
    }
    index = (hashnum % 101) % partitionSize;
    return index;
}

int hashFunc::SumAndPrintItemsInIndex(int index)
{
    int number = NumberOfItemsInIndex(index);


    item* Ptr = HashTable[index];
    if (Ptr->location == "empty")
    {
        //cout<< "index = "<<index<<" is empty";
    }
    else {
        for (int i = 0; i<number; i++)
        {
            if (i ==0)
            HashTable[index]->speed = Ptr->speed; //This stores the total. The tweak here is just so that it doesnt add the first no twice
            else
            HashTable[index]->speed += Ptr->speed; //This stores the total
            Ptr = Ptr->next;
        }

    }
    //cout<<"Sum of elements in SP["<<index<<"] is "<<HashTable[index]->speed<<endl; //This is what we pass to SP i.,e what redies in the smaller holes i.e v
    return HashTable[index]->speed;
    //Habving a 1 D array is better. It makes dividing easier
}
