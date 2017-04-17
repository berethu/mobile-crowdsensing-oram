#pragma once

#include "AES.hpp"

#include <fstream>
#include <sstream>

class File {
public:
	static void Read(std::fstream &file, byte_t *data, size_t len);
	static void Write(std::fstream &file, byte_t *data, size_t len);

	static void ReadNum(std::stringstream &file, byte_t *data, size_t len);


	static size_t GetLength(std::fstream &file);
};
