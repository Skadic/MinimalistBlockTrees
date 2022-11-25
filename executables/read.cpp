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
    uint64_t len;
    ifs.read((char*) &len, sizeof(uint64_t));
    auto* cbt = new CBlockTree(ifs);

    for (size_t i = 0; i < 100; i++) {
      std::cout << (char) cbt->access(i) << std::endl;
    }

    delete cbt;
    return 0;
}

