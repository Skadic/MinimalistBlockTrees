
#include <pointer_based/BlockTree.h>
#include <compressed/CBlockTree.h>

#include <sstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <bitset>
#include <unordered_set>


int main(int argc, char **argv) {

    if (argc < 2) {
      std::cerr << "Please input file" << std::endl;
      exit(1);
    }

    std::string input;
    std::ifstream t(argv[1]);
    std::stringstream buffer;
    buffer << t.rdbuf();
    input = buffer.str();

    std::unordered_set<int> characters;
    for (char c: input) {
        characters.insert(c);
    }

    BlockTree* bt = new BlockTree(input, 2, 32);
    bt->process_back_pointers();
    bt->clean_unnecessary_expansions();
    for (char c: characters)
        bt->add_rank_select_support(c);

    auto* cbt = new CBlockTree(bt);

    std::stringstream ss;
    ss << argv[1];
    ss << ".bt";
    std::string out_path = ss.str();

    uint64_t n = input.length();

    std::ofstream ot(out_path);

    ot.write(reinterpret_cast<const char*>(&n), sizeof(uint64_t));

    cbt->serialize(ot);
    ot.close();

    delete bt;
    delete cbt;
    return 0;
}

