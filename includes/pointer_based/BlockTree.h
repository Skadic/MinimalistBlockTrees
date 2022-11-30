#ifndef BLOCKTREE_PBLOCKTREE_H
#define BLOCKTREE_PBLOCKTREE_H

#include <blocktree.utils/RabinKarp.h>
#include <blocktree.utils/HashString.h>
#include <pointer_based/blocks/BackBlock.h>
#include <pointer_based/blocks/Block.h>

#include <string>
#include <unordered_map>

class BlockTree {
    void block_scan(std::vector<Block*>&, int, std::unordered_map<HashString, std::vector<Block*>>&);
public:
    int r_; // Arity
    int root_arity_;
    int leaf_length_;
    std::string input_; // Input sequence of the Tree
    Block* root_block_;
    bool rank_select_support_;

    BlockTree(std::string& source, int r, int root_block_arity, int leaf_length, bool clean = false, bool rank_select_support = false);
    ~BlockTree();



    int access(int);
    void add_rank_select_support(int);
    int rank(int, int);
    int select(int, int);

    void process_back_pointers_heuristic();
    void process_back_pointers();
    void clean_unnecessary_expansions();

    void process_level_heuristic(std::vector<Block*>&);
    void process_level(std::vector<Block*>&);

    void forward_window_block_scan(std::vector<Block*>& level, int window_size, int N, std::unordered_map<HashString, std::vector<Block*>>& hashtable);
    void forward_pair_window_block_scan(std::vector<Block*>& level, int pair_window_size, int N, std::unordered_map<HashString, std::vector<std::pair<Block*, Block*>>>& pair_hashtable);

    std::vector<Block*> next_level(std::vector<Block*>&);
    // Returns a vector of levels of nodes of the tree where
    // each level is represented by a vector of its nodes (left-to-right).
    //
    // A simple levelwise (left-to-right) traversal of the tree would be:
    //     for (std::vector<Block*> level : bt->levelwise_iterator()) {
    //         for (Block* b : level) {
    //             ...
    std::vector<std::vector<Block*>> levelwise_iterator();


    void print();
};

#endif //BLOCKTREE_PBLOCKTREE_H
