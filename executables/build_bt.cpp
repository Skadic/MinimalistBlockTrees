
#include <pointer_based/BlockTree.h>
#include <compressed/CBlockTree.h>

#include <iostream>
#include <fstream>
#include <sstream>
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

    if (argc < 3) {
      std::cerr << "Please input tree arity (tau)" << std::endl;
      exit(1);
    }

    size_t arity = atoi(argv[2]);

    if (argc < 4) {
      std::cerr << "Please input root arity (s)" << std::endl;
      exit(1);
    }

    size_t root_arity = atoi(argv[3]);

    if (argc < 5) {
      std::cerr << "Please input max leaf length" << std::endl;
      exit(1);
    }

    size_t leaf_length = atoi(argv[4]);

    std::stringstream ss;
    ss << argv[1] << "_arit" << arity << "_root" << root_arity << "_leaf" << leaf_length << ".bt";
    std::string out_path = ss.str();

    std::cout << "building block tree with parameters:" << 
      "\narity: " << arity << 
      "\nroot arity: " << root_arity << 
      "\nmax leaf length: " << leaf_length << 
      "\nsaving to " << out_path << std::endl;

    std::string input;
    std::ifstream t(argv[1]);
    std::stringstream buffer;
    buffer << t.rdbuf();
    input = buffer.str();

    BlockTree* bt = new BlockTree(input, arity, root_arity, leaf_length, true, false);
    auto* cbt = new CBlockTree(bt);

    uint64_t n = input.length();

    std::ofstream ot(out_path);

    ot.write(reinterpret_cast<const char*>(&n), sizeof(uint64_t));

    cbt->serialize(ot);
    ot.close();

    delete bt;
    delete cbt;
    return 0;
}

