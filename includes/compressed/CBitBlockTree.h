#ifndef BLOCKTREE_PCBLOCKTREE_H
#define BLOCKTREE_PCBLOCKTREE_H

#include <pointer_based/BlockTree.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include <unordered_map>


class CBitBlockTree {
public:
    int arity_; // Arity
    int first_level_length_;
    int number_of_levels_;

    std::vector<sdsl::bit_vector*> bt_bv_; // 1 when is Internal Block
    std::vector<sdsl::rank_support_v<1>*> bt_bv_rank_;
    std::vector<sdsl::int_vector<>*> bt_offsets_;
    sdsl::bit_vector* leaf_bv_;


    sdsl::int_vector<>* bt_first_level_prefix_ranks_;

    std::vector<sdsl::int_vector<>*> bt_ranks_;
    std::vector<sdsl::int_vector<>*> bt_second_ranks_;


    CBitBlockTree(BlockTree* source_tree, int one_symbol);
    CBitBlockTree(std::istream& input);
    virtual ~CBitBlockTree();

    int access(int i);
    int rank_0(int i);
    int rank_1(int i);
    int select_0(int rank);
    int select_1(int rank);

    int size();
    int get_partial_size();
    void serialize(std::ostream& output);
};


#endif //BLOCKTREE_PCBLOCKTREE_H
