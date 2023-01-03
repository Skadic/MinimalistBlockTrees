#ifndef BLOCKTREE_PLEAVEBLOCK_H
#define BLOCKTREE_PLEAVEBLOCK_H

#include "Block.h"

class LeafBlock : public Block {
  public:
    LeafBlock(Block *parent, int64_t start_index, int64_t end_index, std::string &source);
    ~LeafBlock();

    int     add_rank_select_support(int character);
    int64_t size();

    int rank(int character, int i);
    int select(int character, int rank);

    int access(int);
};

#endif // BLOCKTREE_PLEAVEBLOCK_H
