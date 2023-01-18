#ifndef BLOCKTREE_PBACKBLOCK_H
#define BLOCKTREE_PBACKBLOCK_H

#include "Block.h"

///
/// @brief A block pointing back at another block as its source.
///
/// This block starts reading from its `first_block` at its given `offset`.
/// It might be that this blocks content spans not just one, but two blocks in its source.
/// If that is the case, then `second_block` is populated with that second block to continue reading from.
///
class BackBlock : public Block {
  public:
    ///
    /// @brief Create a new back block.
    ///
    /// @param parent This block's parent.
    /// @param start_index The start index of this block in the input.
    /// @param end_index The (inclusive) end index of this block in the input.
    /// @param input The input string.
    /// @param first_block The first block from which to copy.
    /// @param second_block The (potentially null) second block from which to copy.
    /// @param offset The offset inside the first block from which to start reading.
    ///
    BackBlock(Block             *parent,
              int64_t            start_index,
              int64_t            end_index,
              const std::string &input,
              Block             *first_block,
              Block             *second_block,
              int                offset);
    ~BackBlock() override;

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
    [[nodiscard]] int select(int character, int i) const override;
};

#endif // BLOCKTREE_PBACKBLOCK_H
