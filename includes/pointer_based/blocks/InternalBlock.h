#ifndef BLOCKTREE_PLAZYINTERNALBLOCK_H
#define BLOCKTREE_PLAZYINTERNALBLOCK_H

#include "Block.h"

class InternalBlock : public Block {
public:

    InternalBlock(Block* parent, int64_t start_index, int64_t end_index, std::string& source);
    ~InternalBlock();

    std::vector<Block*>& children(int leaf_length, int arity);
    void clean_unnecessary_expansions();

    bool is_leaf();
    int access(int i);
    int add_rank_select_support(int character);

    int rank(int character, int i);
    int select(int character, int rank);

};

#endif //BLOCKTREE_PLAZYINTERNALBLOCK_H
