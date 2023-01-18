#ifndef BLOCKTREE_PLEAVEBLOCK_H
#define BLOCKTREE_PLEAVEBLOCK_H

#include "Block.h"

///
/// @brief A leaf block which does not have any children.
///
class LeafBlock : public Block {
  public:
    ///
    /// @brief Create a new leaf block.
    ///
    /// @param parent This block's parent.
    /// @param start_index This block's start index in the input.
    /// @param end_index This block's (inclusive) end index in the input.
    /// @param input The input string.
    ///
    LeafBlock(Block *parent, int64_t start_index, int64_t end_index, const std::string &input);
    ~LeafBlock() override;

    ///
    /// @brief Add support for `rank` and `select` to this block, for the given character.
    ///
    /// @param character A character to add rank select support for.
    /// @return The number of times this character exists in this block.
    ///
    int add_rank_select_support(const int character) override;

    ///
    /// @brief Returns this block's size.
    ///
    /// @return The number of characters this block spans.
    ///
    int64_t size() const;

    ///
    /// @brief A rank query on this block for a given character.
    ///
    /// On a leaf block, this just scans until index i and returns the count.
    ///
    /// @param character A character
    /// @param i An index
    /// @return The number of times the character occurs up to and including the given index.
    ///
    [[nodiscard]] int rank(int character, int i) const override;

    ///
    /// @brief A select query on this block for a given character.
    ///
    /// On a leaf block this just scans until it finds the right character.
    ///
    /// @param character The character to search for
    /// @param rank The rank of the character to look for
    /// @return The position of the rank-th occurrence of the given character
    ///
    [[nodiscard]] int select(int character, int rank) const override;

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
};

#endif // BLOCKTREE_PLEAVEBLOCK_H
