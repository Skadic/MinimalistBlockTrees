

#include <memory>
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

    auto bt = std::unique_ptr<BlockTree>(new BlockTree(input, arity, root_arity, leaf_length, true, true));
    auto cbt = std::unique_ptr<CBlockTree>(new CBlockTree(&*bt));

    uint64_t n = input.length();


    while (true) {
      std::cout << "> ";
      size_t i;
      std::cin >> i;
      
      std::cout << "bt[" << i << "] = " << (char) cbt->access(i) << std::endl;
    }
    
    return 0;
}

