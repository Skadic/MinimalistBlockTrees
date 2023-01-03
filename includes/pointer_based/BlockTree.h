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
    /// @brief Scans through this level's blocks and inserts each block into the hash table.
    ///
    /// The hashes are Rabin-Karp fingerprints of the section of the input that the block spans.
    ///
    /// @param level A list of block pointers that represents a level in the tree
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-Hashed string to all blocks that match its hash.
    ///
    void block_scan(std::vector<Block *> &, int, std::unordered_map<HashString, std::vector<Block *>> &);

  public:
    int         arity_; // Arity
    int         root_arity_;
    int         leaf_length_;
    std::string input_; // Input sequence of the Tree
    Block      *root_block_;
    bool        rank_select_support_;

    /**
     * @sigma_ Alphabet Size
     */
    uint8_t alphabet_size_;

    /**
     * @to_ascii_ Maps an entry in the block tree to its corresponding extended ascii character
     */
    std::array<uint8_t, 256> to_ascii_;

    /**
     * @to_alphabet_ Maps an extended ascii character to its representation in this block tree
     */
    std::array<uint8_t, 256> to_alphabet_;

    BlockTree(std::string &source,
              int          r,
              int          root_block_arity,
              int          leaf_length,
              bool         clean               = false,
              bool         rank_select_support = false);
    ~BlockTree();

    int  access(int);
    void add_rank_select_support(int);
    int  rank(int, int);
    int  select(int, int);

    void process_back_pointers_heuristic();

    ///
    /// @brief Actually build the block tree and insert back-links for repeated blocks.
    ///
    void process_back_pointers();
    void clean_unnecessary_expansions();

    void process_level_heuristic(std::vector<Block *> &);
    void process_level(std::vector<Block *> &);

    ///
    /// @brief Scan for occurrences of strings (of size window_size) that also exist in the hashtable and set the first
    /// occurrences and backlinks of the blocks that match these strings.
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
    void forward_window_block_scan(std::vector<Block *>                                 &level,
                                   int                                                   window_size,
                                   int                                                   N,
                                   std::unordered_map<HashString, std::vector<Block *>> &hashtable);

    ///
    /// @brief Scan for occurrences of strings (of size window_size) that also exist in the hashtable
    /// and for all pairs of blocks that match these strings they are marked.
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
    void forward_pair_window_block_scan(
        std::vector<Block *>                                                     &level,
        int                                                                       pair_window_size,
        int                                                                       N,
        std::unordered_map<HashString, std::vector<std::pair<Block *, Block *>>> &pair_hashtable);

    ///
    /// @brief This takes a vector of blocks and returns a vector of their children in order.
    ///
    /// As this method's name says, this is supposed to be used with input vectors
    /// which consist of blocks of a single level in the tree from left to right.
    ///
    /// @param level A vector of pointers to blocks whose children to return.
    ///
    std::vector<Block *> next_level(std::vector<Block *> &);

    // Returns a vector of levels of nodes of the tree where
    // each level is represented by a vector of its nodes (left-to-right).
    //
    // A simple levelwise (left-to-right) traversal of the tree would be:
    //     for (std::vector<Block*> level : bt->levelwise_iterator()) {
    //         for (Block* b : level) {
    //             ...
    std::vector<std::vector<Block *>> levelwise_iterator();

    void print();
};

#endif // BLOCKTREE_PBLOCKTREE_H
