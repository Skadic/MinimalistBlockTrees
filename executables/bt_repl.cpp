

#include <compressed/CBlockTree.h>
#include <cstdlib>
#include <pointer_based/BlockTree.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <ranges>
#include <sstream>
#include <string>

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

    if (argc < 6) {
        std::cerr << "Please input the size of saves prefixes and suffixes" << std::endl;
        exit(1);
    }
    size_t prefix_suffix_size = atoi(argv[5]);

    std::cout << "building block tree with parameters:"
              << "\narity: " << arity << "\nroot arity: " << root_arity << "\nmax leaf length: " << leaf_length
              << "\nprefix_suffix_size: " << prefix_suffix_size << std::endl;

    std::string       input;
    std::ifstream     t(argv[1]);
    std::stringstream buffer;
    buffer << t.rdbuf();
    input = buffer.str();

    auto bt = std::unique_ptr<BlockTree>(new BlockTree(input, arity, root_arity, leaf_length, true, false));
    bt->add_fast_substring_support(prefix_suffix_size);
    auto cbt = std::unique_ptr<CBlockTree>(new CBlockTree(&*bt));

    while (true) {
        std::cout << "> ";
        std::string input;
        std::cin >> input;
        if (std::ranges::all_of(input.begin(), input.end(), [](char c) { return std::isdigit(c); })) {
            const size_t i = atoi(input.c_str());
            std::cout << "bt[" << i << "] = " << (char) cbt->access(i) << std::endl;
        } else {
            const auto sep = input.find(":");
            if (sep == std::string::npos) {
                std::cout << "<integer> = char random access\n<integer>:<integer> = substring" << std::endl;
                continue;
            }

            const auto l = atoi(input.substr(0, sep).c_str());
            const auto r = atoi(input.substr(sep + 1).c_str());

            std::cout << "l: " << l << ", r: " << r << std::endl;

            if (l > r) {
                std::cout << "left bound cannot be greater than right bound" << std::endl;
                continue;
            }

            auto buf = new char[r - l + 2];
            std::fill(buf, buf + r - l + 2, 0);
            cbt->substr(buf, l, r - l + 1);
            std::cout << "bt[" << l << ":" << r << "] = " << buf << std::endl;
            delete[] buf;
        }
    }

    return 0;
}
