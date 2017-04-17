// sPASS system with HE (U) as the top layer with each containing V ORAM


#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>   
#include <gmp.h>
#include <cstdlib>
#include <system_error>
#include <assert.h>
#include <ios>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <queue>
#include <istream>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <unordered_map>  
#include <time.h>
#include <cstdio>
#include <cstring>
#include <mutex>
std::mutex myMutex;
#include <thread>
#include "general.c"
#include "perm_b64.h"
#include "perm_bas.c"
//HE
extern "C"{
    #include <paillier.h>
}

#include <math.h>
//ORAM
#include "oramNew/src/ORAM.hpp"
#include "oramNew/src/Timer.hpp"
#include "oramNew/src/FileSystem.hpp"
#include "oramNew/src/File.hpp"
#include "oramNew/src/RAMStore.hpp"
#include "oramNew/src/Log.hpp"
#include "oramNew/src/ArraySystem.hpp"
#include "HashFunc/hashFunc.h"
#include "bloom_filter.hpp"

static const int noOfBits = 64;
static const int noOfUsers = 50000;
static const int secMemSize = 512;

 // Optimization and parameter tweaking

bool OSPoptimization(true);
bool dynamicStrategy(true);
bool staticStrategy(false);
bool posMapOptimization(true);
bool memSchedulingOptimization(true);
//Overflowing SP
int OSPval = 30;


t_bool error;

// Type used to store and retrieve codes.
using CodeType = std::uint16_t;

namespace globals {

// Dictionary Maximum Size (when reached, the dictionary will be reset)
const CodeType dms {std::numeric_limits<CodeType>::max()};

} // namespace globals



std::ifstream::pos_type filesize(string filename);
double AggregateRequest(size_t treeDepth, char* input, std::string output, std::string stashbin, std::string posbin);
int MBtoDepth(size_t mb, size_t blockSize = 4096);
double Profile(std::function<void()> fun);
void WriteNumToORAM(FileSystem &files, char* number, std::string outputFilename);
void DynamicCommit(int U, int V, std::vector<char*>& encSum, int& lambdaMaxm);
void StaticCommit(int U, int V, std::vector<char*>& encSum, int N);
void EvictOSP(vector<char*>& matrix, int U, int V, int& lambdaMax);
void set_error();
void odometer(t_longint loop);
void benesNetwork1(const tr_benes* self, const ta_index c_tgt);
// void benesNetwork1a(const tr_benes* self, const ta_index c_tgt);
// void benesNetwork2(const tr_benes* self, const ta_index c_tgt);
// void benesNetwork3(const tr_benes* self);
void bitIndex();
int getVal (int value);
void shuffleBenes();
void benesNetworkPerm();
void benesParity();
std::vector<char> operator + (std::vector<char> vc, char c);
void compress(std::istream &is, std::ostream &os);
void decompress(std::istream &is, std::ostream &os);
void print_usage(const std::string &s = "", bool su = true);


enum StoreType {
    MEM_TYPE,
    ORAM_TYPE
};


using namespace std ;

// from a matrix, retrieve an nxn slice starting at row i, col j
// throws std::out_of_range on error
class partitioner
{

public:

template< typename T >
vector<T> change2Dto1Darray(const vector< vector<T> >& matrix, int numrows, int numcols)
{
    int w = 0;
    int n = numrows * numcols;
vector<T> v(n);
for (int i = 0; i < numrows; i++)
for (int j =0; j<numcols; j++)
v[w++] = matrix[i][j];

return v;
}


template< typename T >
vector<T> partition_slicer( const vector< vector<T> >& matrix, std::size_t m,std::size_t n,
                               std::size_t i, std::size_t j )
{
    //const std::size_t N = matrix.size() ;
    vector< vector<T> > slice(n) ;
    vector<T> b; 

    for( std::size_t k = 0 ; k < m ; ++k )
        for( std::size_t l = 0 ; l < n ; ++l )
            slice[k].push_back( matrix.at(k+i).at(l+j) ) ;
    
    vector <int> v = change2Dto1Darray(slice, m,n );
           
    return v ;
}

template< typename T >
vector<T> partition_slicer1D( const vector<T>& matrix, int u, int start)
{

    // u is the number of partitions we wanna have from a given vector
    // start is the piosition it should start pushing from
    //const std::size_t N = matrix.size() ;
    int part = matrix.size()/u;
    //cout<<"This is "<<part<<endl;;
    vector<T> b;
    for(int i = 0; i<part; i++)
    {
        b.push_back(matrix[start]);
        start++;
        if(start > matrix.size())
            return b;
        
    }            
    return b ;
}



//To cross check
template< typename T >
vector< vector<T> > get_slice( const vector< vector<T> >& matrix, std::size_t n,
                               std::size_t i, std::size_t j )
{
    const std::size_t N = matrix.size() ;
    vector< vector<T> > slice(n) ;
    vector<T> b; 

    for( std::size_t k = 0 ; k < n ; ++k )
        for( std::size_t l = 0 ; l < n ; ++l )
            slice[k].push_back( matrix.at(k+i).at(l+j) ) ;
           

    return slice ;
}



//Printing 2D
template< typename T > void print( const vector< vector<T> >& matrix )
{
    for( std::size_t i = 0 ; i < matrix.size() ; ++i )
    {
        for( std::size_t j = 0 ; j < matrix[i].size() ; ++j )
        {
          std::cout << matrix[i][j] << ' ' ;
            }
        std::cout << '\n' ;
    }
    std::cout << "--------------------\n" ;
}

//printing 1D
template< typename T > int countPrint( const vector<T>& matrix )
{
    std::cout << '\n' ;
    int count = 0;
    for( std::size_t i = 0 ; i < matrix.size() ; ++i )
    {
             //     std::cout << matrix[i] << ' ' ;
              //    std::cout << '\n' ;
                if (matrix[i] != 0)
                    count++;
    }
    //std::cout << "--------------------\n" ;
return count;
}

};

hashFunc Hashy;
partitioner Partitions;
std::vector<char*> allPartition(Hashy.pSize); //holds all data

std::ofstream simulationResult;
std::string simulationRes = "simulationRes"; //Saving simulation result

int rowA = 0, colA = 0;

string arr[44589400][3]; //Gets data from dataset


int main(int argc, char **argv)
{

if (staticStrategy && dynamicStrategy)
{
    print_usage(std::string("Both staticStrategy and dynamicStrategy cannot be turned on"));
        return EXIT_FAILURE;
}
else if (!staticStrategy && !dynamicStrategy)
{
    print_usage(std::string("Both staticStrategy and dynamicStrategy cannot be turned off"));
        return EXIT_FAILURE;
}

std::string dataSource = "trip_resultNew00.txt";
int U = strtol(argv[1], nullptr, 10);; 
int V = allPartition.size() / U;
int N = U * V;
//Holds the SP
std::vector<int>superPartition[U];
//counts amount of request in a partiton
int count[U];
//aggregates the data to partitions
vector<int> aggregator(Hashy.pSize);
//This is the array that will be getting the ciphertext value i.e encSum
std::vector<char*>encSum(Hashy.pSize);
std::vector<string>posMapVector;


simulationRes += std::to_string(U);
simulationRes += ".txt";

simulationResult.open(simulationRes, std::ios_base::app);

//Hashing time for distribution of data
    int rows = Hashy.getData(dataSource, rowA,  colA, arr);
 if (noOfUsers < 100000)
    {
        //More user = more data. The current data is assumed to  be for 100,000 users now
        //Incrementing no of users
        rows = (rows * noOfUsers) / 100000;
    }

    for (int i = 0; i<rows; i++)
    {
      //Even distribution of data into different partitions
        Hashy.AddItem((arr[i][0]+arr[i][1]), stod(arr[i][2]));
    }


   if(memSchedulingOptimization)
    aggregator.reserve(secMemSize);

for(int i = 0; i<secMemSize; i++)//This prints and assigns sum
{
aggregator[i] = Hashy.SumAndPrintItemsInIndex(i);
//cout<<aggregator[i]<<endl; //this is passed to the SP
}


auto start = std::chrono::high_resolution_clock::now();

std::cout<<std::fixed;
 simulationResult<<std::fixed;
 simulationResult<<std::setprecision(7);
    std::cout<<std::setprecision(7);

//encryption time (HE)
auto startHE = std::chrono::high_resolution_clock::now();

paillier_pubkey_t* pub;//The public key
       paillier_prvkey_t* prv;//The private key 

    paillier_keygen(1024, &pub, &prv,paillier_get_rand_devurandom);

    paillier_ciphertext_t* SPencrypt[Hashy.pSize];
    paillier_plaintext_t* SPdecrypt; //For debugging purpose only: decrypts cipher text. I can make this an array but its not needed since decryption is not what we are into in this project

      paillier_plaintext_t* SPplaintext[Hashy.pSize];

      //Getting plain text of all SP element
     for(int i = 0; i<Hashy.pSize; i++)
     {
     SPplaintext[i]=paillier_plaintext_from_ui(aggregator[i]);
     //encryption function
     SPencrypt[i] = paillier_enc(0, pub, SPplaintext[i], paillier_get_rand_devurandom); //encrypted summation
    //This is what will be passed to the aggregate requesst
    }


//Encrypted Sum
for(int i = 0; i<Hashy.pSize; i++)
{
encSum[i] = mpz_get_str(NULL, 10, SPencrypt[i]->c);
 cout<<"Encrypted Sum for partition "<<i <<" "<<encSum[i]<<endl;
}

//End HE
auto endHE = std::chrono::high_resolution_clock::now();

auto tHE = (endHE - startHE).count() * ((double) std::chrono::high_resolution_clock::period::num / std::chrono::high_resolution_clock::period::den);
    // std::cout << "tHE: " <<  tHE << std::endl;
    simulationResult << "tHE: " <<  tHE<<  "  secs"<<std::endl;
    simulationResult << "noOfUsers: "<< noOfUsers<<std::endl;
//cout<<encSum.size()<<endl;;


//Splitting the aggregator to get lambda 
//Spliting time

for(int i = 0, v = 0; i<U; i++, v+=V)
{

superPartition[i] =  Partitions.partition_slicer1D(aggregator,U,v);
count[i] = Partitions.countPrint(superPartition[i]);

}

//maximum lambda for ORAM
int lambdaMax = *std::max_element(count,count+U);
// simulationResult << "lambdaMax : " <<  lambdaMax <<std::endl;

 

auto startORAM = std::chrono::high_resolution_clock::now();

if(OSPoptimization){
EvictOSP(encSum, U,V, lambdaMax);
}
simulationResult << "lambdaMax : " <<  lambdaMax <<std::endl;
if(dynamicStrategy)
{
DynamicCommit(U,V,encSum, lambdaMax);
}
else if (staticStrategy)
{
     StaticCommit(U, V, encSum,N);
}

auto endORAM = std::chrono::high_resolution_clock::now();
//End ORAM

auto end = std::chrono::high_resolution_clock::now();


 auto tORAM = (endORAM - startORAM).count() * ((double) std::chrono::high_resolution_clock::period::num / std::chrono::high_resolution_clock::period::den);

    // std::cout << "tORAM: " <<  tORAM << std::endl;
    simulationResult << "tORAM: " <<  tORAM<< "  secs"<<std::endl;



//Benes network
  if (!init_general()) {
    return 1;
    }

  error = false;
  bitIndex();
  shuffleBenes();
 benesParity(); 
 benesNetworkPerm();


if(posMapOptimization)
{

 enum class Mode {
        Compress,
        Decompress
    };

    Mode m;
    m = Mode::Compress;

    const std::size_t buffer_size {1024 * 1024};

    // these custom buffers should be larger than the default ones
    const std::unique_ptr<char[]> input_buffer(new char[buffer_size]);
    const std::unique_ptr<char[]> output_buffer(new char[buffer_size]);

    std::ifstream input_file;
    std::ofstream output_file;

    input_file.rdbuf()->pubsetbuf(input_buffer.get(), buffer_size);
    // input_file.open(argv[2], std::ios_base::binary);
    input_file.open("posbin.txt", std::ios_base::binary);

    if (!input_file.is_open())
    {
        print_usage(std::string("PosMap `") + argv[2] + "' could not be opened.");
        return EXIT_FAILURE;
    }

    output_file.rdbuf()->pubsetbuf(output_buffer.get(), buffer_size);
    // output_file.open(argv[3], std::ios_base::binary);
    output_file.open("posbin", std::ios_base::binary);


    if (!output_file.is_open())
    {
        print_usage(std::string("posMap `") + argv[3] + "' could not be opened.");
        return EXIT_FAILURE;
    }

    try
    {
        input_file.exceptions(std::ios_base::badbit);
        output_file.exceptions(std::ios_base::badbit | std::ios_base::failbit);

        if (m == Mode::Compress)
            compress(input_file, output_file);
        else
        if (m == Mode::Decompress)
            decompress(input_file, output_file);
    }
    catch (const std::ios_base::failure &f)
    {
        print_usage(std::string("PosMap input/output failure: ") + f.what() + '.', false);
        return EXIT_FAILURE;
    }
    catch (const std::exception &e)
    {
        print_usage(std::string("Caught exception: ") + e.what() + '.', false);
        return EXIT_FAILURE;
    }

}
else{
    std::ifstream positionMapData("posbin.txt");
std::stringstream posMapBuffer;
posMapBuffer << positionMapData.rdbuf();
string posMapData = posMapBuffer.str();
std::cout<<"\n";
std::cout<<"Writing position Map to memory";
for(int i = 0; i<posMapData.length(); i++)
{
posMapVector.push_back(posMapData);
}

}

auto finalTime = (end - start).count() * ((double) std::chrono::high_resolution_clock::period::num / std::chrono::high_resolution_clock::period::den);
    std::cout << "Experiment Completed" << std::endl;
    std::cout << "Total Time taken: " <<  finalTime << std::endl;
    simulationResult << "Total Time taken: " <<  finalTime<< std::endl;

auto fileSize = filesize(dataSource);

simulationResult << "Transfered data size: "<<fileSize<< " bytes"<< std::endl;

auto bandwidth = fileSize / finalTime;

simulationResult << "Bandwidth is "<<bandwidth<<" bytes per sec"<< std::endl;

return 0;
}






double AggregateRequest(size_t treeDepth, char* input, std::string output, std::string stashbin, std::string posbin)
{
    //std::lock_guard<std::mutex> guard(myMutex);
    //All for ORAM
    Timer timer;
    timer.Start();
    auto start = std::chrono::high_resolution_clock::now();
    AES::Setup();

    srand(time(NULL));
    
    bytes<Key> key {"RandomKeyIshSomeLongStringToEnc"}; //32bytes string for encryption
    
    // Retrive the depth of the tree
    size_t depth = treeDepth;//strtol(argv[1], nullptr, 10);
    printf("depth = %zu\n", depth);

    size_t blockSize = 1024 * (noOfBits/8);
    size_t blockCount = Z*(pow(2, depth + 1) - 1);
    
    printf("block count = %zu, block size = %zu\n", blockCount, blockSize);

    size_t storeBlockSize = 0;
    size_t storeBlockCount = 0;

    // What's backing the secure store
    StoreType storeType = MEM_TYPE;
    
    // Parse the secure store type (ORAM)
    StoreType secureType = ORAM_TYPE;

    switch (secureType) {
        case ORAM_TYPE:
            storeBlockSize = IV + AES::GetCiphertextLength(Z*(sizeof (int32_t) + blockSize));
            storeBlockCount = blockCount/4;
            break;
        default:
            Log::Write(Log::FATAL, "Secure store type missing");
    }

    // Create the store
    BlockStore *store = nullptr;
    
    switch (storeType) {
        case MEM_TYPE:
            store = new RAMStore(storeBlockCount, storeBlockSize);
            break;

        default:
            Log::Write(Log::FATAL, "Store type missing");
    }
    
    // Create the secure store
    BlockStore *secureStore = nullptr;

    switch (secureType) {
        case ORAM_TYPE:
            secureStore = new ORAM(store, depth, blockSize, key, stashbin, posbin);
            break;

        default:
            Log::Write(Log::FATAL, "Secure store type missing");
    }

    printf("#blocks = %zu\n", secureStore->GetBlockCount());

    
    // if (!secureStore->WasSerialised()) {
    //     puts("zeroing blocks");

    //     for (size_t i = 0; i < secureStore->GetBlockCount(); i++) {
    //         printf("\rintialising %zu/%zu", i+1, secureStore->GetBlockCount());
    //         fflush(stdout);

    //         block b(blockSize, i % 256);
    //         secureStore->Write(i, b);
    //     }
    //     puts("");
    // }
    
    // files.Load();
    FileSystem files(secureStore);


    //Start iteration here

    files.LoadNum();
        
    Profile([&]() {
            if (!files.AddNum(input)) {
            WriteNumToORAM(files, input, output);
        }
    });
    
   files.SaveNum();
    

    //block b(secureStore->GetBlockSize(), 17);
    //secureStore->Write(0, b);

    //block b = secureStore->Read(0);
    //printf("block contains %d\n", (int) b[0]);

    
    
    assert(secureStore != NULL);
delete secureStore;
secureStore = NULL;
    assert(store != NULL);
delete store;
store = NULL;

    AES::Cleanup();
    timer.Stop();
    auto end = std::chrono::high_resolution_clock::now();
    double elapsedTime = timer.GetElapsedTime();
    
    printf("Elapsed time = %f\n", elapsedTime);
    std::cout << "Time taken: " << (end - start).count() * ((double) std::chrono::high_resolution_clock::period::num / std::chrono::high_resolution_clock::period::den) << std::endl;

}



void WriteNumToORAM(FileSystem &files, char* number, std::string outputFilename)
{
    BlockStore *store = files.GetBlockStore();

    // Retrieve metadata
    FileInfo info = files.GetFileInfo(number); 

    // Open output file
    std::fstream file;
    file.open(outputFilename, std::ios::out | std::ios::binary | std::ios::trunc);

    size_t blockSize = store->GetBlockSize();

    for (size_t i = 0; i < info.length; i += blockSize) {
        size_t writeLength = std::min(blockSize, info.length - i);
        size_t pos = info.blocks[i/blockSize];

        // Read from the store
        block buffer = store->Read(pos);  //This gets the file from oram store or position map

        // Write it to file
        file.write((char *) buffer.data(), writeLength);

        printf("\r%zu / %zu", i/blockSize + 1, info.length/blockSize);
        fflush(stdout);
    }
    puts("\n");

    file.close();
}



int MBtoDepth(size_t mb, size_t blockSize /*= 4096*/)
{
    mb *= 1024*1024;
    
    size_t blocks = ceil(mb/(double) blockSize);
    size_t buckets = ceil(blocks/(double) Z);
    
    return ceil(log2(buckets+1))-2;
}

double Profile(std::function<void()> fun)
{
    Timer timer;
    timer.Start();
    
    fun();
    
    timer.Stop();
    double elapsedTime = timer.GetElapsedTime();
    
    //printf("Elapsed time = %f\n", elapsedTime);
        
    return elapsedTime;
}



struct lookaheadData {
std::string CurrentReq;    
double TamtReq;
    
    struct lookaheadData *nextReq;
};


void DynamicCommit(int U, int V, std::vector<char*>& encSum, int& lambdaMaxm){

vector<char*>encryptedValue[U];
size_t depth;
int KLookahead = 100;
int lambdaReq = 0;int lambdaMax = 0;
double TamtReq = 0, TamtInitial = 0;
std::string reqInitial;
std::vector<lookaheadData> lookaheads;
lookaheadData lookahead;

   bloom_parameters parameters;

   // estimated no of data to be inserted
   parameters.projected_element_count = 10000000;
   parameters.false_positive_probability = 0.0001;
   parameters.random_seed = 0xA5A5A5A5;

   parameters.compute_optimal_parameters();
   bloom_filter filter(parameters);

for(int i = 0, v = 0; i<U; i++, v+=V)
{
encryptedValue[i] =  Partitions.partition_slicer1D(encSum,U,v);
}

simulationResult<<"U: "<<U<<endl;
// simulationResult<<"lamdaMax "<<lambdaMaxm<<endl; 

depth = MBtoDepth(noOfBits);//depth of tree
 

int totalCount = 0;  
int optN = lambdaMaxm * U;  
// simulationResult<<"optN "<<optN <<endl;
// simulationResult<<"\n ";

         for(int SPu = 0; SPu <U; SPu++)
          for (int SPv = 0; SPv <lambdaMaxm; SPv++) {

if (filter.contains(encryptedValue[SPu][SPv]) == 0)
    {
        //add the data into bloomFilter
    filter.insert(encryptedValue[SPu][SPv]);
    lambdaReq++;
    }

 if (lambdaReq > lambdaMax)
{
lambdaMax = lambdaReq;
if (TamtReq < TamtInitial)
{
    reqInitial = encryptedValue[SPu][SPv]; 
        lookaheads.clear();
}
else{

		lookahead.CurrentReq = encryptedValue[SPu][SPv];
		lookahead.TamtReq = TamtReq; 

        lookaheads.push_back (lookahead);
        if(lookaheads.size() == KLookahead)
        {
            if (!lookaheads.empty())
            reqInitial = lookaheads.back().CurrentReq;
        }
}
}


totalCount++;
std::string output = "output";
//output += std::to_string(SPu); 
output += ".txt";

std::string stashbin = "stashbin";
//stashbin += std::to_string(SPu);
stashbin += ".txt";

std::string posbin = "posbin";
//posbin += std::to_string(SPu);
posbin += ".txt";

// #pragma omp barrier
cout<<totalCount<<" out of "<<optN<<endl;
AggregateRequest(depth,encryptedValue[SPu][SPv], output, stashbin, posbin);
}

}



void  StaticCommit(int U, int V, std::vector<char*>& encSum, int N){
vector<char*>encryptedValue[U];
size_t depth;
int count[U];

for(int i = 0, v = 0; i<U; i++, v+=V)
{

encryptedValue[i] =  Partitions.partition_slicer1D(encSum,U,v);
count[i] = Partitions.countPrint(encryptedValue[i]);

}
int lambdaMaxEnc = *std::max_element(count,count+U); 

 simulationResult<<"lambdaMaxEnc "<<lambdaMaxEnc<<endl; 
simulationResult<<"U "<<U<<endl;

    depth = 13; //gotten based on experiments
 

int totalCount = 0; 
simulationResult<<"N"<<N <<endl;

         for(int SPu = 0; SPu <U; SPu++)
          for (int SPv = 0; SPv <V; SPv++) {


totalCount++;
std::string output = "output";
//output += std::to_string(SPu); 
output += ".txt";

std::string stashbin = "stashbin";
//stashbin += std::to_string(SPu);
stashbin += ".txt";

std::string posbin = "posbin";
//posbin += std::to_string(SPu);
posbin += ".txt";

// #pragma omp barrier
cout<<totalCount<<" out of "<<N<<endl;
AggregateRequest(depth,encryptedValue[SPu][SPv], output, stashbin, posbin);
}

}



void EvictOSP(vector<char*>& matrix, int U, int V, int& lambdaMax) //partition = encripted value matrix = encSum
{
std::vector<char*> partition[U];
std::priority_queue<char*> priorityQH;
int count[U];
int noOfOSP;
int lambdaEvict;
double x;
srand (time(NULL));
    int lambdaCount = 0;
    if(OSPval > 20)
    {
        OSPval = OSPval % 20;
        // simulationResult << "OSPval : " <<  OSPval <<std::endl;
        if (OSPval > 20)
            OSPval = rand() % 8 + 13;
        if (OSPval <1)
            OSPval = 1;
         // simulationResult << "OSPvalFinal : " <<  OSPval <<std::endl;
    }

          for(int i = 0, v = 0; i<OSPval; i++, v+=V)
{

partition[i] =  Partitions.partition_slicer1D(matrix,U,v);
count[i] = Partitions.countPrint(partition[i]);

}
 lambdaEvict = *std::max_element(count,count+U);
lambdaEvict = getVal(lambdaEvict);

if (lambdaEvict < lambdaMax){
lambdaMax = lambdaEvict;
}

   for( std::size_t i = 0 ; i < U ; ++i )
    {
matrix.pop_back();
                if (matrix[i] != 0)
                {
                    priorityQH.push(matrix[i]);
                    lambdaCount++;
                    
                        while (!priorityQH.empty())
                {
                    lambdaEvict--;
                    priorityQH.pop();
                }
                    if(lambdaCount > lambdaMax){
                         lambdaCount = lambdaEvict;
                         } 
                }

    }
}


void set_error()
{
  error = true;
  // exit(1);
  }

 void odometer(t_longint loop) {
  if ((loop & 0x1fff) == 0) {
    printf("%i    \x0d",loop);
    fflush(0);
    }
  }


 void benesNetwork1(const tr_benes* self, const ta_index c_tgt) {

  t_int q;
  t_bits mask;
  t_bits x;
  ta_index c_inv_tgt;

  invert_perm(c_tgt,c_inv_tgt);

  mask = 1;
  for (q=0; q<=benesBits-1; ++q) {
    x = benes_fwd(self,mask);
    if (x != lo_bit << c_inv_tgt[q]) {
      printf("Error: benesNetwork1\n");
      set_error();
      }
    mask = mask << 1;
    }
  }

 // void benesNetwork1a(const tr_benes* self, const ta_index c_tgt) {

 //  t_int q;
 //  t_bits mask;
 //  t_bits x;

 //  mask=1;
 //  for (q=0; q<=benesBits-1; ++q) {
 //    x=benes_fwd(self,lo_bit << c_tgt[q]);
 //    if (x != mask) {
 //      printf("Error: benesNetwork1a\n");
 //      set_error();
 //      }
 //    mask = mask << 1;
 //    }
 //  }

 // void benesNetwork2(const tr_benes* self, const ta_index c_tgt) {

 //  t_int q;
 //  t_bits mask;
 //  t_bits x;

 //  mask = 1;
 //  for (q=0; q<=benesBits-1; ++q) {
 //    x=benes_bwd(self,mask);
 //    if (x != lo_bit << c_tgt[q]) {
 //      printf("Error: benesNetwork2\n");
 //      set_error();
 //      }
 //    mask = mask << 1;
 //    }
 //  }

 // void benesNetwork3(const tr_benes* self) {

 //  t_int q;
 //  t_bits mask;
 //  t_bits x;

 //  mask = 1;
 //  for (q=0; q<=benesBits-1; ++q) {
 //    x = benes_fwd(self,benes_bwd(self,mask));
 //    if (x != mask) {
 //      printf("Error: benesNetwork3\n");
 //      set_error();
 //      }
 //    mask = mask << 1;
 //    }
 //  }


 void bitIndex() {

  t_longint loop;
  t_subword j,k;
  t_bits x,y;

  printf("Bit index\n");
  for (loop = noOfUsers; loop > 0; --loop) {
    odometer(loop);
    x = random_bits();
    y = x;
    for (j = 0; j<=ld_bits-1; ++j) {
      x = bit_index_complement(x,j);
      }
    if (x != general_reverse_bits(y,benesBits-1)) {
      printf("bit_index_complement: Error\n");
      set_error();
      }

    y = x;
    for (j = ld_bits-1; j>=1; --j) {
      x = bit_index_swap(x,j,j-1);
      }
    if (x != shuffle(y,0,ld_bits)) {
      printf("bit_index_swap: Error\n");
      set_error();
      }

    y = x;
    j = random_int(ld_bits);  // 0..ld_bits-1
    k = random_int(ld_bits);  // 0..ld_bits-1
    x = bit_index_swap_complement(x,j,k);
    if (x != bit_index_complement(
            bit_index_complement(
              bit_index_swap(y,j,k),
            j),
          k)) {
      printf("bit swap: Error\n");
      set_error();
      }
    }
  }

 void shuffleBenes() {

  t_longint loop;
  t_int i;
  t_bits x,y,z;
  t_subword sw_entities,ld_col,ld_row;

  for (loop = noOfUsers; loop > 0; --loop) {
    odometer(loop);
    x = random_bits();
    do {
      ld_row = random_int(ld_bits);  // 0..ld_bits-1
      ld_col = random_int(ld_bits);  // 0..ld_bits-1
      sw_entities = random_int(ld_bits);  // 0..ld_bits-1
      } while (!(sw_entities+ld_row+ld_col <= ld_bits));
    y = x;
    for (i = 1; i<=(t_int)(ld_col); ++i) {
      y = shuffle(y, sw_entities,ld_row+ld_col+sw_entities);
      }
    z = transpose(x, sw_entities,ld_row,ld_col);
    if (y != z) {
      printf("transpose: Error\n");
      set_error();
      }
    z = shuffle_power(x,sw_entities,ld_row+ld_col+sw_entities,ld_col);
    if (y != z) {
      printf("shuffle_power: Error\n");
      set_error();
      }
    z = unshuffle_power(z,sw_entities,ld_row+ld_col+sw_entities,ld_col);
    if (x != z) {
      printf("unshuffle_power: Error\n");
      set_error();
      }
    }
  }


 void benesNetworkPerm() {

  t_longint loop;
  tr_benes benes;
  ta_index c_tgt;
  t_int i;


  for (loop = noOfUsers; loop > 0; --loop) {
    odometer(loop);
    random_perm(c_tgt);
    gen_benes(&benes,c_tgt);
    if (benes.b1.cfg[0] != 0) {
      printf("BN: Error (stage 0)\n");
      set_error();
      }
    for (i=0; i<=ld_bits-1; ++i) {
      if ((benes.b1.cfg[i] &
          ~(a_bfly_mask[i] & ~a_bfly_lo[i+1])) != 0) {
        printf("BN: Error (Waksman)\n");
        set_error();
        }
      }
    benesNetwork1(&benes,c_tgt);
    // benesNetwork1a(&benes,c_tgt);
    // benesNetwork2(&benes,c_tgt);
    // benesNetwork3(&benes);
    }
    printf("Permutations complete\n");
  }

int getVal (int value)
{
    std::ostringstream os;
os << value;
std::string valStr = os.str();
std::string val1 = valStr.substr (0,2);
std::string val2 = valStr.substr (2,1);
int finalVal = stoi(val1) + stoi(val2);
return finalVal;
}

 void benesParity() {

  t_longint loop;
  tr_benes benes;
  ta_index c_tgt;
  t_bool p1, p2;
  t_int i,j;
  t_int x;

  // printf("Parity of a Benes network\n");
  identity_perm(c_tgt);
  p1 = false;
  for (loop = noOfUsers; loop > 0; --loop) {
    odometer(loop);
    gen_benes(&benes,c_tgt);  // need target indexes for each c_tgt
    p2 = benes_parity(&benes);
    if (p1 != p2) {
      printf("benesParity: Error\n");
      set_error();
      }
    // exchange two different entries
    do {
      i = random_int(benesBits);  // 0..benesBits-1
      j = random_int(benesBits);  // 0..benesBits-1
      } while (!(i != j));
    x = c_tgt[i];
    c_tgt[i] = c_tgt[j];
    c_tgt[j] = x;
    p1 = ! p1;
    }
  }


//  Helper hash functor, to be used on character containers.
struct HashCharVector {

    // returns Hash of the vector of characters.
    std::size_t operator () (const std::vector<char> &vc) const
    {
        return std::hash<std::string>()(std::string(vc.cbegin(), vc.cend()));
    }
};


// returns vector resulting from appending `c` to `vc`
std::vector<char> operator + (std::vector<char> vc, char c)
{
    vc.push_back(c);
    return vc;
}

//  Compresses the contents of `is` and writes the result to `os`.

void compress(std::istream &is, std::ostream &os)
{
    std::unordered_map<std::vector<char>, CodeType, HashCharVector> dictionary;

    // reset the dictionary to its initial contents
    const auto reset_dictionary = [&dictionary] {
        dictionary.clear();

        const long int minc = std::numeric_limits<char>::min();
        const long int maxc = std::numeric_limits<char>::max();

        for (long int c = minc; c <= maxc; ++c)
        {
            // to prevent Undefined Behavior, resulting from reading and modifying
            // the dictionary object at the same time
            const CodeType dictionary_size = dictionary.size();

            dictionary[{static_cast<char> (c)}] = dictionary_size;
        }
    };

    reset_dictionary();

    std::vector<char> s;
    char c;

    while (is.get(c))
    {
        // dictionary's maximum size was reached
        if (dictionary.size() == globals::dms)
            reset_dictionary();

        s.push_back(c);

        if (dictionary.count(s) == 0)
        {
            // to prevent Undefined Behavior, resulting from reading and modifying
            // the dictionary object at the same time
            const CodeType dictionary_size = dictionary.size();

            dictionary[s] = dictionary_size;
            s.pop_back();
            os.write(reinterpret_cast<const char *> (&dictionary.at(s)), sizeof (CodeType));
            s = {c};
        }
    }

    if (!s.empty())
        os.write(reinterpret_cast<const char *> (&dictionary.at(s)), sizeof (CodeType));
}

//  Decompresses the contents of `is` and writes the result to `os`.
void decompress(std::istream &is, std::ostream &os)
{
    std::vector<std::vector<char>> dictionary;

    // reset the dictionary to its initial contents
    const auto reset_dictionary = [&dictionary] {
        dictionary.clear();
        dictionary.reserve(globals::dms);

        const long int minc = std::numeric_limits<char>::min();
        const long int maxc = std::numeric_limits<char>::max();

        for (long int c = minc; c <= maxc; ++c)
            dictionary.push_back({static_cast<char> (c)});
    };

    reset_dictionary();

    std::vector<char> s;
    CodeType k; // Key

    while (is.read(reinterpret_cast<char *> (&k), sizeof (CodeType)))
    {
        // dictionary's maximum size was reached
        if (dictionary.size() == globals::dms)
            reset_dictionary();

        if (k > dictionary.size())
            throw std::runtime_error("invalid compressed code");

        if (k == dictionary.size())
            dictionary.push_back(s + s.front());
        else
        if (!s.empty())
            dictionary.push_back(s + dictionary.at(k).front());

        os.write(&dictionary.at(k).front(), dictionary.at(k).size());
        s = dictionary.at(k);
    }

    if (!is.eof() || is.gcount() != 0)
        throw std::runtime_error("corrupted compressed file");
}


//  Prints usage information and a custom error message.
void print_usage(const std::string &s, bool su)
{
    if (!s.empty())
        std::cerr << "\nERROR: " << s << '\n';

    if (su)
    {
        std::cerr << "\nInvalid usage\n";
    }

    std::cerr << std::endl;
}

//get filesize of data sent btw ORAM and HE
std::ifstream::pos_type filesize(string filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg(); 
}