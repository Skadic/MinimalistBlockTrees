#ifndef BLOCKTREE_PBLOCK_H
#define BLOCKTREE_PBLOCK_H

#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

///
/// @brief A block inside the block tree spanning a given span of the input text.
///
class Block {
  public:
    /// The parent block
    Block *parent_;

    /// The index in the input text at which this block starts
    int64_t start_index_;

    /// The index in the input text at which this block starts
    int64_t end_index_;

    /// The string this block tree is built from
    const std::string &input_;

    /// Maps a character to the number of times it appears in this block
    std::unordered_map<int, int> pop_counts_;

    /// This is a `BackBlock`-specific field. `BackBlock`s point to a source to the left of them from which they copy
    /// characters. This source which may span one or two blocks. This stores how many times character c appears in the
    /// part of the source that lies in the *first block* it copies from.
    std::unordered_map<int, int> pop_counts_in_first_block_;

    /// The first block of this block's source (if this is a back block)
    Block *first_block_;

    /// The second block of this block's source (if this is a back block).
    /// This is only populated if the source is at an offset where the source does not fit into only first_block
    Block *second_block_;

    /// The offset at which this block's source starts in its first_block
    int offset_{};

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

    /// Needed for fast substring support. This stores the prefix of the string represented by this block
    /// represented.
    /// The size is given by the `prefix_suffix_size_` field in `BlockTree`.
    /// This field is truncated if it would extend past the end of the input string!
    std::string_view prefix_;

    /// Needed for fast substring support. This stores the prefix of the string represented by this block
    /// represented.
    /// The size is given by the `prefix_suffix_size_` field in `BlockTree`.
    /// This field is truncated if it would extend past the end of the input string!
    std::string_view suffix_;

    ///
    /// @brief Create a new block.
    ///
    /// @param parent This block's parent.
    /// @param start_index Ths block's start index in the input.
    /// @param end_index Ths block's start index in the input.
    /// @param input The input string.
    ///
    Block(Block *parent, int64_t start_index, int64_t end_index, const std::string &input);
    virtual ~Block();

    ///
    /// @brief This block's length.
    ///
    /// @return The number of characters this block spans.
    ///
    int64_t length() const;

    ///
    /// @brief Returns a copy of the segment of the string which this block represents.
    ///
    std::string represented_string() const;

    ///
    /// @brief Add support for `rank` and `select` to this block, for the given character.
    ///
    /// @param character A character to add rank select support for.
    /// @return The number of times this character exists in this block.
    ///
    virtual int add_rank_select_support(int character);

    ///
    /// @brief Add support for fast substring queries.
    ///
    /// This populates the `prefix_suffix_` field.
    ///
    virtual void add_fast_substring_support(int prefix_suffix_size);

    ///
    /// @brief Accesses the given index.
    ///
    /// @param i An index
    /// @return The character at index i in the input string.
    ///
    virtual int access(int i) const;

    ///
    /// @brief Get a substring from the input text.
    ///
    /// @param buf The char buffer to write to
    /// @param index The index to start reading from.
    /// @param len The length of the string to read.
    /// @return The pointer position after writing
    ///
    virtual char *substr(char *buf, int index, int len) const;

    ///
    /// @brief A rank query on this block for a given character.
    ///
    /// Note, that this requires rank support on this block.
    ///
    /// @param character A character
    /// @param i An index
    /// @return The number of times the character occurs up to and including the given index.
    ///
    virtual int rank(int character, int i) const;

    ///
    /// @brief A select query on this block for a given character.
    ///
    /// Note, that this requires select support on this block.
    ///
    /// @param character The character to search for
    /// @param rank The rank of the character to look for
    /// @return The position of the rank-th occurrence of the given character
    ///
    virtual int select(int character, int rank) const;

    ///
    /// @brief Returns this block's children.
    ///
    /// @param leaf_length The tree's leaf length
    /// @param arity This block's arity
    ///
    virtual std::vector<Block *> &children(int leaf_length, int arity);

    ///
    /// @brief If an internal block has only leaves as its children, no other blocks are pointing to it and this block
    /// itself points back to another block, then replace it with a back block pointing back at its source.
    ///
    virtual void clean_unnecessary_expansions();

    ///
    /// @brief Replaces one of this block's children with another block.
    ///
    /// @param old_child The old child to be replaced.
    /// @param new_child The new child to take its place
    ///
    void replace_child(Block *old_child, Block *new_child);

    ///
    /// @brief Checks whether this block is a leaf block.
    ///
    /// @return true, if this block is a leaf block *or* a back block, false if it is an internal block.
    ///
    virtual bool is_leaf() const;

    ///
    /// @brief Get the first `count` chars from this Block's prefix.
    ///
    /// This method returns at most as many characters as are actually saved in this block.
    ///
    /// @param count The number of characters.
    ///
    std::string_view prefix(size_t count) const;

    ///
    /// @brief Get the last `count` chars from this Block's suffix.
    ///
    /// This method returns at most as many characters as are actually saved in this block.
    ///
    /// @param count The number of characters.
    ///
    std::string_view suffix(size_t count) const;
};

#endif // BLOCKTREE_PBLOCK_H
