#ifndef BLOCKTREE_PBACKBLOCK_H
#define BLOCKTREE_PBACKBLOCK_H

#include "Block.h"

class BackBlock : public Block {
  public:
    BackBlock(Block       *parent,
              int64_t      start_index,
              int64_t      end_index,
              std::string &source,
              Block       *first_block,
              Block       *second_block,
              int          offset);
    ~BackBlock();

    int access(int i);
    int add_rank_select_support(int character);

    int rank(int character, int i);
    int select(int character, int i);
};

#endif // BLOCKTREE_PBACKBLOCK_H
