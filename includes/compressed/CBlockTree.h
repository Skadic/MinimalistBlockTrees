#ifndef BLOCKTREE_PCBLOCKTREE_H
#define BLOCKTREE_PCBLOCKTREE_H

#include <pointer_based/BlockTree.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include <unordered_map>
#include <unordered_set>

class CBlockTree {
  public:
    /// This arity of each block in this tree, except for the root node.
    int arity_;
    /// The root block's arity.
    int root_arity_;
    /// The block size of the blocks in the first level of this tree.
    int first_level_block_size_;
    /// The number of levels in this tree.
    int number_of_levels_;
    /// Whether this tree is built with rank select support.
    bool rank_select_support_;
    /// The number of characters in the input string;
    size_t input_size_;
    /// The size of the prefixes and suffixes of a block's represented string. This is saved for every block larger than
    /// this value, allowing fast substring queries.
    size_t prefix_suffix_size_;

    /// For each level l, stores a bitvector B,
    /// where B[j] == 1 if the block at index j on level l is an internal block
    /// This maps level -> block index -> is internal
    std::vector<sdsl::bit_vector *> is_internal_;

    /// For each level, this stores a rank data structure for the `is_internal` bit vector
    /// This maps level -> block index -> rank
    std::vector<sdsl::rank_support_v<1> *> is_internal_ranks_;

    /// For each level, stores an intvector which for each non-internal block contains the offset from the start of the
    /// string to its source.
    /// For example:
    /// ABA|CAC|BAC|CBA
    ///
    /// For the block representing "CBA", this would store 5, since it is copying from index 5.
    ///
    /// This maps level -> block index (skipping internal blocks) -> offset
    std::vector<sdsl::int_vector<> *> offsets_;

    /// The string represented by the lowest level of this tree.
    /// This is transformed by the `mapping_`
    sdsl::int_vector<> *leaf_string_;

    /// Contains the first and last prefix_suffix_size_ symbols of each node of each level.
    /// One int_vector contains the first prefix_suffix_size_ symbols of node 0,
    /// then the last prefix_suffix_size_ symbols of node 0, then the first log_sigma(n) symbols of node 1, and so on.
    std::vector<sdsl::int_vector<> *> prefix_suffix_symbols_;

    /// Maps an integer code number to the character it represents
    sdsl::int_vector<> *alphabet_;
    /// Maps a character to an integer code number
    std::unordered_map<char, int> mapping_;

    /// For the first level (corresponds to the lowest level in the original tree with no missing blocks) of the tree,
    /// saves the ranks up to (and not including) a specific block for each character.
    /// This maps character -> block index -> rank
    std::unordered_map<int, sdsl::int_vector<> *> first_level_ranks_;

    /// This stores the number of times a character appears inside of each of the tree's blocks.
    /// This maps character -> level -> block index -> pop count
    std::unordered_map<int, std::vector<sdsl::int_vector<> *>> pop_counts_;

    /// This is specific to back blocks and stores the number of times a character appears in the part of this
    /// block's source which is contained in the first block it spans.
    /// For example:
    /// ABAA|CADC|AACA|BACA
    ///
    /// The block "AACA" copies from index 2 to (inclusively) index 5. the first "AA" is contained in the first block of
    /// its source. So it would store 2 for character A, since it appears two times in the part of its source that is in
    /// the first block.
    ///
    /// This maps character -> level -> block index (skipping internal blocks) -> pop count
    std::unordered_map<int, std::vector<sdsl::int_vector<> *>> first_block_pop_counts_;

    explicit CBlockTree(BlockTree *source_tree);
    explicit CBlockTree(std::istream &input);
    virtual ~CBlockTree();

    ///
    /// @brief Get the character at the given index in the input.
    ///
    /// @param index An index in the input text.
    /// @return The character at the given index in the input text.
    ///
    [[nodiscard]] int access(int i) const;

    /// @brief Copy the substring at an index and of the given length into a buffer.
    ///
    /// @param buf The buffer to write to.
    /// @param index The index to read from.
    /// @param len The length of the substring to copy
    /// @return The position of the buffer after the written string
    char *substr(char *buf, size_t index, size_t len) const;

    ///
    /// @brief Gets the number the given character's occurrences up to an index.
    ///
    /// @param character A character to count.
    /// @param index The index up to which occurences are counted.
    /// @return The number of times the character appears up to this index.
    ///
    [[nodiscard]] int rank(int character, int index) const;

    ///
    /// @brief Gets the index of the rank-th occurrence of the given character.
    ///
    /// @param character A character to search for.
    /// @param rank The desired rank of the character.
    /// @return An index i, such that `access(i) == c` and `rank(c, i) == rank`.
    ///
    [[nodiscard]] int select(int character, int rank) const;

    ///
    /// @brief Gets the approximate in-ram size of this data structure in RAM.
    ///
    /// @return The size in bytes.
    ///
    [[nodiscard]] int size() const;
    [[nodiscard]] int get_partial_size() const;
    void              serialize(std::ostream &character) const;

  private:
    [[nodiscard("for internal calls the buffer should be adjusted")]] char *
    substr_internal(char *buf, size_t index, size_t len) const;
};

#endif // BLOCKTREE_PCBLOCKTREE_H
