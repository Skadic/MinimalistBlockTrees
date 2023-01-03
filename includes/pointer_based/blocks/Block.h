#ifndef BLOCKTREE_PBLOCK_H
#define BLOCKTREE_PBLOCK_H

#include <string>
#include <unordered_map>
#include <vector>

class Block {
  public:
    /// The parent block
    Block *parent_;

    /// The index in the input text at which this block starts
    int64_t start_index_;

    /// The index in the input text at which this block starts
    int64_t end_index_;

    /// The string this block tree is built from
    std::string &source_;

    /// Maps a character to the number of times it appears in this block
    std::unordered_map<int, int> ranks_;

    /// TODO Idek
    std::unordered_map<int, int> second_ranks_;

    /// The first block of this block's source (if this is a back block)
    Block *first_block_;
    /// The second block of this block's source (if this is a back block).
    /// This is only populated if the source is at an offset where the source does not fit into only first_block
    Block *second_block_;
    /// The offset at which this block's source starts in its first_block
    int offset_;

    /// This will be set to true if the string represented by this block and its right neighbor
    /// was found somewhere else in the text
    bool left_;

    /// This will be set to true if the string represented by this block and its left neighbor
    /// was found somewhere else in the text
    bool right_;

    /// The number of blocks that are pointing to this block,
    /// i.e. the number of blocks whose source is (at least partly) contained in this block.
    int pointing_to_me_;

    /// The index this block has in its respective level
    int level_index_;

    /// The level index of first_block_
    int first_occurrence_level_index_;

    /// Pointers to this block's children
    std::vector<Block *> children_;

    Block(Block *parent, int64_t start_index, int64_t end_index, std::string &source);
    virtual ~Block();

    ///
    /// @brief This block's length.
    ///
    /// @return The number of characters this block spans.
    ///
    int64_t     length() const;
    std::string represented_string();

    virtual int add_rank_select_support(int character);

    virtual int access(int i);
    virtual int rank(int character, int i);
    virtual int select(int character, int rank);

    virtual std::vector<Block *> &children(int leaf_length, int arity);
    virtual void                  clean_unnecessary_expansions();
    void                          replace_child(Block *old_child, Block *new_child);

    ///
    /// @brief Checks whether this block is a leaf block.
    ///
    /// @return true, if this block is a leaf block *or* a back block, false if it is an internal block.
    ///
    virtual bool is_leaf();
};

#endif // BLOCKTREE_PBLOCK_H
