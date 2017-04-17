#include "ArraySystem.hpp"
#include "BlockStore.hpp"
#include "File.hpp"
#include "Log.hpp"

#include <sstream>
#include <algorithm>

ArraySystem::ArraySystem(BlockStore *store)
: store(store), files()
{}

ArraySystem::~ArraySystem()
{}

int ArraySystem::GetAvailableID()
{
	static size_t id = 0;
	
	if (id + 1 >= store->GetBlockCount()) {
		Log::Write(Log::FATAL, "out of space");
	}
	
	return id++;
}

bool ArraySystem::Add(double arrayname[], int size)
{
	// Keep track of file info
	std::string arraynam; //I dont know why i did this
	ArrayInfo &info = files[arraynam];
	info.length = size; 
	
	size_t blockSize = store->GetBlockSize();
		
		for (size_t i = 0; i < info.length; i += blockSize)
		{

		block b(blockSize, 0);
		//real array implementation

			b.push_back(arrayname[i]); //copies the value of the array into the vector
		
		// Generate random blockID
		size_t bid = GetAvailableID();
		info.blocks.push_back(bid);
	
		store->Write(bid, b);	
		
		//printf("\r%zu / %zu", i/CHUNK + 1, info.length/CHUNK);
		//fflush(stdout);
	}
	//puts("");
	
	
	return true;
}

bool ArraySystem::Remove(std::string filename) //no need
{
	if (files.find(filename) == files.end()) {
		Log::Write(Log::WARNING, "file not in database");
		return false;
	}
	
	// TODO: This doesn't actually free the IDs

	files.erase(filename);
	
	return true;
}

std::pair<std::string, ArrayInfo> ArraySystem::LoadArrayInfo(std::string line)
{		
	std::istringstream sstream(line);

	std::string filename;
	sstream >> filename;
	
	ArrayInfo info;
	sstream >> info.length;

	size_t blockSize = store->GetBlockSize();

	for (size_t i = 0; i < info.length; i += blockSize) {
		size_t bid; //This is the bucket ID
		sstream >> bid;
	
		info.blocks.push_back(bid);
	}
	
	return std::make_pair(filename, info);
}

void ArraySystem::Load()
{
	std::ifstream file("filemapArray.txt");
	
	if (!file) {
		return;
	}
	
	std::string line;

	while (std::getline(file, line)) {
	
		files.insert(LoadArrayInfo(line));
	}
	
	file.close();
}

std::string ArraySystem::SaveArrayInfo(std::string filename, ArrayInfo info)
{
	std::ostringstream sstream;
	
	sstream << filename;
	sstream << ' ';
	sstream << info.length;

	size_t blockSize = store->GetBlockSize();
	
	for (size_t i = 0; i < info.length; i += blockSize) {
		sstream << ' ';
		sstream << info.blocks[i/blockSize];
	}
	
	return sstream.str();
}

void ArraySystem::Save()
{
	std::ofstream file("filemapArray.txt");
	
	for (auto f : files) {
		file << SaveArrayInfo(f.first, f.second);
		file << '\n';
	}
	
	file.close();
}

BlockStore *ArraySystem::GetBlockStore()
{
	return store;
}

ArrayInfo ArraySystem::GetArrayInfo(std::string filename)
{
	return files[filename];
}
