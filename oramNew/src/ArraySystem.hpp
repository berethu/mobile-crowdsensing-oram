#pragma once

#include <cstddef>
#include <unordered_map>
#include <list>
#include <vector>

struct ArrayInfo {
	size_t length;
	std::vector<size_t> blocks;
};

class BlockStore;

class ArraySystem {
	BlockStore *store; //We can make this shared pointyer to allow for dynamism
	
	std::unordered_map<std::string, ArrayInfo> files;
	
	int GetAvailableID();

	std::pair<std::string, ArrayInfo> LoadArrayInfo(std::string line);
	std::string SaveArrayInfo(std::string filename, ArrayInfo info);

public:
	ArraySystem(BlockStore *store);
	~ArraySystem();

	bool Add(double arrayname[], int size);
	bool Remove(std::string filename);
	
	void Load();
	void Save();
	
	BlockStore *GetBlockStore();

	ArrayInfo GetArrayInfo(std::string filename);
};
