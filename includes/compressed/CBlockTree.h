#ifndef BLOCKTREE_PCBLOCKTREE_H
#define BLOCKTREE_PCBLOCKTREE_H

#include <pointer_based/BlockTree.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include <unordered_map>

class CBlockTree {
  public:
    /// This arity of each block in this tree, except for the root node.
    int arity_;
    /// The root block's arity.
    int root_arity_;
    /// The block size of the blocks in the lowest complete level of this tree.
    int lowest_complete_level_block_size_;
    /// The number of levels in this tree.
    int number_of_levels_;
    /// Whether this tree is built with rank select support.
    bool rank_select_support_;

    /// For each level l, stores a bitvector B,
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

    /// @brief For the lowest complete level (i.e. with no blocks missing) of the tree, saves the ranks up to (and not
    /// including) a specific block for each character This maps character -> block index -> rank
    std::unordered_map<int, sdsl::int_vector<> *> lowest_complete_level_ranks_;

    /// @brief This stores the number of times a character appears inside of each of the tree's blocks.
    /// This maps character -> level -> block index -> pop count
    std::unordered_map<int, std::vector<sdsl::int_vector<> *>> pop_counts_;
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
