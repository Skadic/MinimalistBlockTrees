#include <pointer_based/BlockTree.h>
#include <compressed/CBlockTree.h>

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>


int main(int argc, char **argv) {

    if (argc < 2) {
      std::cerr << "Please input file" << std::endl;
      exit(1);
    }

    if (!std::filesystem::exists(argv[1])) {
      std::cerr << "File " << argv[1] << " does not exist" << std::endl;
      exit(1);
    }

    std::ifstream ifs(argv[1]);
    auto* cbt = new CBlockTree(ifs);

    std::ofstream ofs(std::string(argv[1]) + ".dec");
    auto buf = new char[cbt->input_size_];
    cbt->substr(buf, 0, cbt->input_size_);
    ofs << buf;
    
    delete[] buf;
    delete cbt;
    return 0;
}

