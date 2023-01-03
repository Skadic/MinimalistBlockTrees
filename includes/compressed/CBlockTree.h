#ifndef BLOCKTREE_PCBLOCKTREE_H
#define BLOCKTREE_PCBLOCKTREE_H

#include <pointer_based/BlockTree.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include <unordered_map>

class CBlockTree {
  public:
    int  arity_; // Arity
    int  root_arity_;
    int  lowest_level_block_length_;
    int  number_of_levels_;
    bool rank_select_support_;

    /// @brief For each level l, stores a bitvector B, 
    /// where B[j] == 1 if the block at index j on level l is an internal block 
    /// This maps level -> block index -> is internal
    std::vector<sdsl::bit_vector *>        is_internal_;
    std::vector<sdsl::rank_support_v<1> *> bt_bv_rank_;
    std::vector<sdsl::int_vector<> *>      bt_offsets_;
    sdsl::int_vector<>                    *leaf_string_;

    /// @brief Contains the first and last log_sigma(n) symbols of each node of each level
    /// one int_vector contains the first log_sigma(n) symbols of node 0, then the last log_sigma(n) symbols of node 0,
    /// then the first log_sigma(n) symbols of node 1, and so on.
    std::vector<sdsl::int_vector<> *> start_end_symbols_;

    sdsl::int_vector<>           *alphabet_;
    std::unordered_map<char, int> mapping_;

    /// @brief For the lowest level of the tree, saves the ranks up to (and not including) a specific block for each character
    /// This maps character -> block index -> rank  
    std::unordered_map<int, sdsl::int_vector<> *> lowest_complete_level_prefix_ranks_;

    /// @brief For the lowest level this stores the ranks of this character inside the current block
    /// This maps character -> block index -> rank 
    std::unordered_map<int, std::vector<sdsl::int_vector<> *>> block_ranks_;
    std::unordered_map<int, std::vector<sdsl::int_vector<> *>> bt_second_ranks_;

    CBlockTree(BlockTree *source_tree);
    CBlockTree(std::istream &input);
    virtual ~CBlockTree();

    int access(int i);
    int rank(int character, int i);
    int select(int character, int rank);

    int  size();
    int  get_partial_size();
    void serialize(std::ostream &character);
};

#endif // BLOCKTREE_PCBLOCKTREE_H
