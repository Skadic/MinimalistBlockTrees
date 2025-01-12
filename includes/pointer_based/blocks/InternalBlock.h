#ifndef BLOCKTREE_PLAZYINTERNALBLOCK_H
#define BLOCKTREE_PLAZYINTERNALBLOCK_H

#include "Block.h"

///
/// @brief An internal block in the block tree. This block neither points back at another block, nor is it a leaf with
/// no children.
///
class InternalBlock : public Block {
  public:
    ///
    /// @brief Create a new internal block.
    ///
    /// @param parent This block's parent.
    /// @param start_index Ths block's start index in the input.
    /// @param end_index Ths block's start index in the input.
    /// @param input The input string.
    ///
    InternalBlock(Block *parent, int64_t start_index, int64_t end_index, const std::string &input);
    ~InternalBlock() override;

    ///
    /// @brief Returns this block's children.
    ///
    /// @param leaf_length The tree's leaf length.
    /// @param arity This block's arity
    ///
    std::vector<Block *> &children(int leaf_length, int arity) override;

    ///
    /// @brief If an internal block has only leaves as its children, no other blocks are pointing to it and this block
    /// itself points back to another block, then replace it with a back block pointing back at its source.
    ///
    void clean_unnecessary_expansions() override;

    ///
    /// @brief Checks whether this block is a leaf. Since this is an internal block, this is false.
    ///
    /// @return false
    ///
    bool is_leaf() const override;

    ///
    /// @brief Accesses the given index.
    ///
    /// @param i An index
    /// @return The character at index i in the input string.
    ///
    [[nodiscard]] int access(int i) const override;

    ///
    /// @brief Get a substring from the input text.
    ///
    /// @param buf The char buffer to write to
    /// @param index The index to start reading from.
    /// @param len The length of the string to read.
    /// @return The pointer position after writing
    ///
    char *substr(char *buf, int index, int len) const override;

    ///
    /// @brief Add support for `rank` and `select` to this block, for the given character.
    ///
    /// @param character A character to add rank select support for.
    /// @return The number of times this character exists in this block.
    ///
    int add_rank_select_support(int character) override;

    ///
    /// @brief A rank query on this block for a given character.
    ///
    /// Note, that this requires rank support on this block.
    ///
    /// @param character A character
    /// @param i An index
    /// @return The number of times the character occurs up to and including the given index.
    ///
    [[nodiscard]] int rank(int character, int i) const override;

    ///
    /// @brief A select query on this block for a given character.
    ///
    /// Note, that this requires select support on this block.
    ///
    /// @param character The character to search for
    /// @param rank The rank of the character to look for
    /// @return The position of the rank-th occurrence of the given character
    ///
    [[nodiscard]] int select(int character, int rank) const override;
};

#endif // BLOCKTREE_PLAZYINTERNALBLOCK_H
