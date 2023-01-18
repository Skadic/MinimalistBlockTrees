#ifndef BLOCKTREE_PBLOCKTREE_H
#define BLOCKTREE_PBLOCKTREE_H

#include <array>
#include <blocktree.utils/HashString.h>
#include <blocktree.utils/RabinKarp.h>
#include <pointer_based/blocks/BackBlock.h>
#include <pointer_based/blocks/Block.h>

#include <string>
#include <unordered_map>

class BlockTree {

  public:
    using Level        = std::vector<Block *>;
    using BlockMap     = std::unordered_map<HashString, std::vector<Block *>>;
    using BlockPairMap = std::unordered_map<HashString, std::vector<std::pair<Block *, Block *>>>;

  private:
    ///
    /// @brief Scans through this level block by block. For each block, this hashes its represented string and inserts
    /// it into the hashtable, mapping from the hashed string to all pairs of blocks that have the same hash.
    ///
    /// The hashes are Rabin-Karp fingerprints of the section of the input that the block spans.
    ///
    /// @param level A list of block pointers that represents a level in the tree
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-hashed string to all blocks that match its hash.
    ///
    void hash_blocks(const Level &level, const int N, BlockMap &hashtable) const;

    ///
    /// @brief Scans through this level block by block with a window of two blocks. hashes the represented string for
    /// each block pair and inserts it into a hashtable, mapping from the hashed string to all pairs of blocks that have
    /// the same hash.
    ///
    /// The hashes are Rabin-Karp fingerprints of the section of the input that the block pair spans.
    ///
    /// @param level A list of block pointers that represents a level in the tree
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-hashed string to all block pairs that match its hash.
    ///
    void hash_block_pairs(const Level &level, const int N, BlockPairMap &hashtable) const;

  public:
    /// This tree's arity (except for the root node)
    const int arity_;
    /// The root node's arity
    const int root_arity_;
    /// The maximum length of a leaf node.
    const int leaf_length_;
    /// The size of the prefixes and suffixes of a block's represented string. This is saved for every internal block,
    /// allowing fast substring queries.
    int prefix_suffix_size_;
    /// Input sequence of the Tree
    const std::string input_;
    /// The root block
    Block *root_block_;

    ///
    /// @brief Create a new block tree from a string with the given properties.
    ///
    /// @param source The source string.
    /// @param arity The arity each node should have except for the root.
    /// @param root_block_arity The root node's arity.
    /// @param leaf_length The maximum length of a leaf node. This is the node size under which recursion stops while
    /// constructing the tree.
    /// @param process_block_tree Whether to process the block tree immediately. Otherwise this can be manually done by
    /// calling `process_block_tree` and `clean_unnecessary_expansions` after.
    /// @param rank_select_support Whether to build this block tree with rank/select support. This can be manually done
    /// by calling `add_rank_select_support` for every desired character.
    ///
    BlockTree(const std::string &source,
              const int          arity,
              const int          root_block_arity,
              const int          leaf_length,
              const bool         process_block_tree  = false,
              const bool         rank_select_support = false);
    ~BlockTree();

    ///
    /// @brief Get the character at the given index in the source.
    ///
    /// @param index An index in the source text.
    /// @return The character at the given index in the source text.
    ///
    int access(const int index) const;

    ///
    /// @brief Get a substring from the input text.
    ///
    /// @param buf The char buffer to write to
    /// @param index The index to start reading from.
    /// @param len The length of the string to read.
    /// @return The pointer position after writing
    ///
    char *substr(char *buf, const int index, const int len) const;

    ///
    /// @brief Add support for fast substring queries.
    ///
    /// This populates the `prefix_suffix_` field in the internal blocks.
    ///
    virtual void add_fast_substring_support(int prefix_suffix_size);

    ///
    /// @brief Adds support for rank and select queries to this block tree.
    ///
    /// @param character The character for which to add the rank and select support.
    ///
    void add_rank_select_support(const int character);

    ///
    /// @brief Gets the number the given character's occurrences up to an index.
    ///
    /// @param character A character to count.
    /// @param index The index up to which occurences are counted.
    /// @return The number of times the character appears up to this index.
    ///
    int rank(const int character, const int index) const;

    ///
    /// @brief Gets the index of the rank-th occurrence of the given character.
    ///
    /// @param character A character to search for.
    /// @param rank The desired rank of the character.
    /// @return An index i, such that `access(i) == c` and `rank(c, i) == rank`.
    ///
    int select(const int character, const int rank) const;

    void process_back_pointers_heuristic();

    ///
    /// @brief Actually build the block tree and insert back-links for repeated blocks.
    ///
    void process_block_tree();

    ///
    /// @brief For each block in the tree, replace unnecessary internal blocks with back blocks where possible. See also
    /// the `InternalBlock` documentation.
    ///
    void clean_unnecessary_expansions();

  private:
    void process_level_heuristic(Level &level);

    ///
    /// @brief Calculate the necessary data, populate the fields in the blocks and create backlinks in this level.
    ///
    /// @param level The level to process.
    ///
    void process_level(Level &level);

    ///
    /// @brief Scan for occurrences of strings (of size window_size) that also exist in the hashtable and set the first
    /// occurrences and backlinks of the blocks that match these strings.
    ///
    /// The hashtable should contain a mapping for each block that exists in the tree. Each mapping should map from the
    /// string represented by the block, to a vector of all blocks that share the same hash value.
    /// This can be calculated with `hash_blocks`.
    ///
    /// This scans through all blocks on this level with a window of the given size and upon recognizing a string that
    /// exists in the hashtable, gets all blocks that match that string and for each block creates a backlink to its
    /// source block(s).
    ///
    /// @param level A list of block pointers that represents a level in the tree.
    /// @param window_size The window size used for scanning.
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-Hashed string to all blocks that match its hash.
    ///
    void window_block_scan(const Level &level, const int window_size, const int N, BlockMap &hashtable) const;

    ///
    /// @brief Scan for occurrences of strings (of size window_size) that also exist in the hashtable
    /// and for all pairs of blocks that match these strings, they are marked.
    ///
    /// The hashtable should contain mappings for each contiguous block pair that exists in the tree.
    /// Each mapping should map from the string represented by the block pair, to a vector of all block pairs that share
    /// the same hash value.
    /// This can be calculated with `hash_block_pairs`.
    ///
    /// This scans through all blocks on this level with a window of the given size and upon recognizing a string that
    /// exists in the hashtable, gets all block pairs that match that string and for each of these pairs marks in both
    /// blocks that the string represented by them both appeared somewhere in the text.
    ///
    /// @param level A list of block pointers that represents a level in the tree.
    /// @param window_size The window size used for scanning.
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-Hashed string to all pairs of blocks that match its hash.
    ///
    void window_block_pair_scan(const Level  &level,
                                const int     pair_window_size,
                                const int     N,
                                BlockPairMap &pair_hashtable) const;

    ///
    /// @brief Create back blocks from the information saved in the blocks on the given level.
    ///
    /// This requires the data in the blocks to be populated. See `window_block_scan` and `window_block_pair_scan` for
    /// more.
    ///
    /// @param level A list of block pointers that represents a level in the tree.
    ///
    void create_back_blocks(Level &level);

  public:
    ///
    /// @brief This takes a vector of blocks and returns a vector of their children in order.
    ///
    /// As this method's name says, this is supposed to be used with input vectors
    /// which consist of blocks of a single level in the tree from left to right.
    ///
    /// @param level A vector of pointers to blocks whose children to return.
    ///
    Level next_level(const Level &level) const;

    /// @brief Returns a vector of levels of blocks of the tree where each level is represented by a vector of its
    /// blocks (left-to-right).
    ///
    /// A simple levelwise (left-to-right) traversal of the tree would be:
    ///     for (std::vector<Block*> level : bt->levelwise_iterator()) {
    ///         for (Block* b : level) {
    ///             ...
    std::vector<Level> levelwise_iterator() const;

    ///
    /// @brief Prints this blocktree to the console
    ///
    void print();
};

#endif // BLOCKTREE_PBLOCKTREE_H
