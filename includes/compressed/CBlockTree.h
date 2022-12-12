#ifndef BLOCKTREE_PCBLOCKTREE_H
#define BLOCKTREE_PCBLOCKTREE_H

#include <pointer_based/BlockTree.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include <unordered_map>


class CBlockTree {
public:
    int arity_; // Arity
    int root_arity_;
    int first_level_length_;
    int number_of_levels_;
    bool rank_select_support_;

    std::vector<sdsl::bit_vector*> bt_bv_; // 1 when is Internal Block
    std::vector<sdsl::rank_support_v<1>*> bt_bv_rank_;
    std::vector<sdsl::int_vector<>*> bt_offsets_;
    sdsl::int_vector<>* leaf_string_;


    sdsl::int_vector<>* alphabet_;
    std::unordered_map<char,int> mapping_;


    std::unordered_map<int,sdsl::int_vector<>*> bt_first_level_prefix_ranks_;

    std::unordered_map<int,std::vector<sdsl::int_vector<>*>> bt_ranks_;
    std::unordered_map<int,std::vector<sdsl::int_vector<>*>> bt_second_ranks_;


    CBlockTree(BlockTree* source_tree);
    CBlockTree(std::istream& input);
    virtual ~CBlockTree();

    int access(int i);
    int rank(int character, int i);
    int select(int character, int rank);

    int size();
    int get_partial_size();
    void serialize(std::ostream& character);
};


#endif //BLOCKTREE_PCBLOCKTREE_H
