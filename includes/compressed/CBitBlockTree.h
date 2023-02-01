#ifndef BLOCKTREE_PCBLOCKTREE_H
#define BLOCKTREE_PCBLOCKTREE_H

#include <pointer_based/BlockTree.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>

#include <unordered_map>

class CBitBlockTree {
  public:
    /// This arity of each block in this tree, except for the root node.
    int arity_;
    /// The block size of the blocks in the first level of this tree.
    int first_level_block_size_;
    /// The number of levels in this tree.
    int number_of_levels_;
    /// The number of bits in the input string;
    size_t input_size_;

    /// For each level l, stores a bitvector B,
    /// where B[j] == 1 if the block at index j on level l is an internal block
    /// This maps level -> block index -> is internal
    std::vector<sdsl::bit_vector *> is_internal_;

    /// For each level, this stores a rank data structure for the `is_internal` bit vector
    /// This maps level -> block index -> rank
    std::vector<sdsl::rank_support_v<1> *> is_internal_ranks_;

    /// For each level, stores an intvector which for each non-internal block contains the offset from the start of the
    /// string to its source.
    /// This maps level -> block index (skipping internal blocks) -> offset
    std::vector<sdsl::int_vector<> *> offsets_;
    /// The block size of the blocks in the first level of this tree.
    // The bit vector resulting in the concatenation of all leaves
    sdsl::bit_vector *leaf_bv_;

    /// For the first level (corresponds to the lowest level in the original tree with no missing blocks) of the tree,
    /// saves the ranks up to (and not including) a specific block.
    /// This maps block index -> rank
    sdsl::int_vector<> *first_level_ranks_;

    /// This stores the number of times a 1 appears inside of each of the tree's blocks.
    /// This maps level -> block index -> pop count
    std::vector<sdsl::int_vector<> *> pop_counts_;

    /// This is specific to back blocks and stores the number of times a character appears in the part of this
    /// block's source which is contained in the first block it spans.
    /// This maps level -> block index (skipping internal blocks) -> pop count
    std::vector<sdsl::int_vector<> *> first_block_pop_counts_;

  
    ///
    /// @brief Creates a new compressed block tree, interpreting one specific character as a 1 and all others as zero.
    ///
    /// @param source_tree The uncompressed block tree to build the compressed tree from.
    /// @param one_symbol The character that should be interpreted as a one.
    ///
    CBitBlockTree(BlockTree *source_tree, int one_symbol);

    ///
    /// @brief Read a compressed block tree from an input stream.
    ///
    /// @param input An input stream supplying a serialized bit block tree.
    ///
    CBitBlockTree(std::istream &input);
    virtual ~CBitBlockTree();

    ///
    /// @brief Get the bit at the given index in the input.
    ///
    /// @param index An index in the input text.
    /// @return The character at the given index in the input.
    ///
    [[nodiscard]] int access(int i) const;

    ///
    /// @brief Gets the number of zeroes occurring up to the given index.
    ///
    /// @param index The index up to which occurences are counted.
    /// @return The number of times 0 appears up to this index.
    ///
    [[nodiscard]] int rank_0(int i) const;

    ///
    /// @brief Gets the number of ones occurring up to the given index.
    ///
    /// @param index The index up to which occurences are counted.
    /// @return The number of times i appears up to this index.
    ///
    [[nodiscard]] int rank_1(int i) const;

    ///
    /// @brief Gets the index of the rank-th occurrence of zero.
    ///
    /// @param rank The desired rank of the zero.
    /// @return An index i, such that `access(i) == 0` and `rank_0(i) == rank`.
    ///
    [[nodiscard]] int select_0(int rank) const;

    ///
    /// @brief Gets the index of the rank-th occurrence of one.
    ///
    /// @param rank The desired rank of the zero.
    /// @return An index i, such that `access(i) == 1` and `rank_1(i) == rank`.
    ///
    [[nodiscard]] int select_1(int rank) const;

    ///
    /// @brief Gets the approximate in-ram size of this data structure in RAM.
    ///
    /// @return The size in bytes.
    ///
    [[nodiscard]] int size() const;
    [[nodiscard]] int get_partial_size() const;
    ///
    /// @brief Serialize the block tree.
    ///
    /// @param output The output stream to which this block tree is written.
    ///
    void              serialize(std::ostream &output) const;
};

#endif // BLOCKTREE_PCBLOCKTREE_H
