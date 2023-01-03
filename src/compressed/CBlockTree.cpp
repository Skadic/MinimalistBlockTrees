#include <compressed/CBlockTree.h>
#include <unordered_set>
#include <vector>


///
/// @brief Gets the lowest level of the given block tree that is not fully comprised of internal blocks.
///
/// In other words, this returns the lowest level of the block tree that is still *complete* in the sense that there are no gaps in the level
/// where blocks are missing.
///
/// @param bt The uncompressed pointer-based block tree
/// @return A vector of pointers pointing to the blocks of the lowest complete level of the given block tree in order from left to right.
///
auto get_lowest_complete_level(BlockTree *bt) -> std::vector<Block*> {
    std::vector<Block*> current_level = {bt->root_block_};
    while (true) {
        for (Block* b: current_level) {
            if(b->is_leaf()) { 
               return current_level;
            }
        }
        current_level = bt->next_level(current_level);
    }

}

///
/// @brief Populates the rank and prefix rank attributes for the lowest complete level of the compressed block tree.
///
/// This corresponds to the fields CBlockTree::lowest_level_prefix_ranks_ and the entry for the lowest level in CBlockTree::lowest_level_block_ranks_.
///
/// @param cbt The compressed block tree whose fields to populate.
/// @param lowest_level_blocks The blocks of the lowest complete level in the original block tree in order from left to right.
///
void populate_lowest_complete_level_ranks(CBlockTree *cbt, const std::vector<Block *> &lowest_complete_level_blocks) {
    /// Saves the current cumulative ranks for each character up to the block in the current iteration
    std::unordered_map<int, int> current_prefix_ranks;

    // Iterate through the rank values for the first block and create a new empty ranks vector for each of them
    // The new rank vectors are used for the compressed block tree
    for (auto &[character, _] : lowest_complete_level_blocks[0]->ranks_) {
        cbt->lowest_complete_level_prefix_ranks_[character] = new sdsl::int_vector<>(lowest_complete_level_blocks.size());
        cbt->block_ranks_[character].push_back(new sdsl::int_vector<>(lowest_complete_level_blocks.size()));
    }

    // iterate through every block of the lowest level and for each character, 
    // calculate the ranks up to this block, as well as the ranks inside this block  
    for (int block_index = 0; block_index < lowest_complete_level_blocks.size(); ++block_index) {
        // For every character, save its prefix rank value up to this block into the new map
        // This will take the prefix rank values calculated in the previous iteration, 
        // since we don't want to save the rank values *before* this block 
        for (auto &[character, rank] : current_prefix_ranks) { 
            (*cbt->lowest_complete_level_prefix_ranks_[character])[block_index] = rank;
        }
      
        // Lookup the ranks for each character inside the current block
        for (auto &[character, _] : lowest_complete_level_blocks[block_index]->ranks_) {
            // Set the ranks of this character for the current block
            (*cbt->block_ranks_[character][0])[block_index] = lowest_complete_level_blocks[block_index]->ranks_[character];
            // Store the cumulative ranks before this block into the temporary map
            // This will be the prefix ranks value for the *next* block (since we only want to cound the ranks *before*) each block
            current_prefix_ranks[character] += lowest_complete_level_blocks[block_index]->ranks_[character];
        }
    }

    // Compress the rank bitvectors
    for (auto &[_, ranks] : cbt->lowest_complete_level_prefix_ranks_) {
        sdsl::util::bit_compress(*ranks);
    }
    for (auto &[_, ranks] : cbt->block_ranks_) {
        sdsl::util::bit_compress(*(ranks[0]));
    }
}

///
/// @brief Populates the rank values inside each block for the given level.
///
/// This modifies the block_ranks_ field. The rank information in the given level 
/// is pushed on top of block_ranks_[character] for each character.
///
/// @param cbt The compressed block tree whose fields to populate.
/// @param level The blocks from the original block tree of the level that is to be processed.
///
void populate_level_block_ranks(CBlockTree *cbt, std::vector<Block *> &level) {
    // Allocate a new int vector for each character to store its ranks
    // We are pushing the new int vector onto the vector containing all intvectors
    // That means that from now on the last element (.back()) for each character is the one we need to worry about for this level
    for (auto &[character, _] : level[0]->ranks_) {
        cbt->block_ranks_[character].push_back(new sdsl::int_vector<>(level.size()));
    }

    // Populate the ranks map with the ranks of the characters inside each block of the current level
    for (int i = 0; i < level.size(); ++i) {
        for (auto &[character, rank] : level[i]->ranks_) {
            (*cbt->block_ranks_[character].back())[i] = rank;
        }
    }

    // Compress the rank bit vectors
    for (auto &[character, ranks] : cbt->block_ranks_) {
        sdsl::util::bit_compress(*ranks.back());
    }
}

CBlockTree::CBlockTree(BlockTree * bt) : arity_(bt->arity_), root_arity_(bt->root_arity_), rank_select_support_(bt->rank_select_support_) {
    std::vector<Block*> lowest_complete_level = get_lowest_complete_level(bt);
    
    /*
    std::cout << "level:" << std::endl;
    for (Block *r : lowest_complete_level) {
      std::cout << "(" << r->start_index_ << "," << r->end_index_ << "), ";    
    }
    std::cout << std::endl;

    std::cout << "first_blocks:" << std::endl;
    for (Block *r : lowest_complete_level) {
      std::cout << "(" << r->first_block_->start_index_ << "," << r->first_block_->end_index_ << "), ";    
    }
    std::cout << std::endl;

    std::cout << "second_blocks:" << std::endl;
    for (Block *r : lowest_complete_level) {
      auto *second_block = r->second_block_;
      if (second_block == nullptr) {
        std::cout << "empty second block" << std::endl;
        break;
      }
      std::cout << "(" << second_block->start_index_ << "," << second_block->end_index_ << "), ";    
    }
    std::cout << std::endl;
    */

    //std::cout << "Second Ranks for 'a':" << std::endl;
    
    /*for (Block *r : lowest_complete_level) {
      std::cout << r->second_ranks_.at('a') << ",";    
    }
    std::cout << std::endl;*/

    // Populate the rank values in the lowest level of the tree that is still complete (i.e. there are no blocks missing)
    populate_lowest_complete_level_ranks(this, lowest_complete_level);

    lowest_level_block_length_ = lowest_complete_level[0]->length();
    number_of_levels_    = 0;

    // Continue further down the tree.
    std::vector<Block *> current_level = lowest_complete_level;
    std::vector<Block *> next_level    = bt->next_level(lowest_complete_level);

    while (next_level.size() != 0) {
        // A bit vector for this level that marks whether the specific block is internal (= 1) or not (= 0).
        sdsl::bit_vector *is_internal_block = new sdsl::bit_vector(current_level.size(), 0);
        
        int number_of_leaves = 0;
        // Iterate through this level and update each block's level index (this block's index inside the current level), 
        // check whether they are internal blocks and count the number of leaves
        for (int i = 0; i < current_level.size(); ++i) {
            current_level[i]->level_index_ = i;

            if (current_level[i]->is_leaf()) {
                (*is_internal_block)[i] = 0;
                ++number_of_leaves;
            } else {
                (*is_internal_block)[i] = 1;
            }
        }

        // Now we populate the rank values in this level
        populate_level_block_ranks(this, next_level);
        
        // TODO Find out what these are
        sdsl::int_vector<> *current_level_offsets = new sdsl::int_vector<>(number_of_leaves);
        std::unordered_map<int, sdsl::int_vector<> *> current_level_second_ranks;

        for (auto [character, _] : current_level[0]->ranks_) {
            current_level_second_ranks[character] = new sdsl::int_vector<>(number_of_leaves);
        }

        // The index of the current block in this level, not counting internal nodes
        int j = 0;
        const int current_length   = current_level.front()->length();
        for (int i = 0; i < current_level.size(); ++i) {
            // If the current block is not internal
            if (!(*is_internal_block)[i]) {
                // This is now the j-th non-internal block on this level
                // Go through the block's "second_rank" for each character 
                // and populate the pre-allocated int vectors with it
                for (auto &[character, second_rank] : current_level[i]->second_ranks_) {
                    (*current_level_second_ranks[character])[j] = second_rank;
                }
                
                (*current_level_offsets)[j++] =
                    current_level[i]->first_block_->level_index_ * current_length + current_level[i]->offset_;
            }
        }

        sdsl::util::bit_compress(*current_level_offsets);
        bt_offsets_.push_back(current_level_offsets);


        for (auto pair : current_level_second_ranks) {
            sdsl::util::bit_compress(*(pair.second));
            (bt_second_ranks_[pair.first]).push_back(pair.second);
        }

        is_internal_.push_back(is_internal_block);
        sdsl::rank_support_v<1> *current_level_bv_rank = new sdsl::rank_support_v<1>(is_internal_block);
        bt_bv_rank_.push_back(current_level_bv_rank);

        current_level = next_level;
        next_level    = bt->next_level(current_level);
        ++number_of_levels_;
    }

    ++number_of_levels_;

    std::vector<Block *> last_level = current_level;

    std::string leaf_string = "";
    for (Block *b : last_level) {
        leaf_string += b->represented_string();
    }

    std::unordered_set<char> alphabet;
    for (char c : leaf_string) {
        alphabet.insert(c);
    }
    alphabet_   = new sdsl::int_vector<>(alphabet.size());
    int counter = 0;
    for (char c : alphabet) {
        mapping_[c]             = counter;
        (*alphabet_)[counter++] = c;
    }
    sdsl::util::bit_compress(*alphabet_);

    leaf_string_ = new sdsl::int_vector<>(leaf_string.size());

    for (int i = 0; i < (*leaf_string_).size(); ++i) {
        (*leaf_string_)[i] = mapping_[leaf_string[i]];
    }

    sdsl::util::bit_compress(*leaf_string_);
}

CBlockTree::CBlockTree(std::istream &in) {
    in.read((char *) &arity_, sizeof(int));
    in.read((char *) &root_arity_, sizeof(int));
    in.read((char *) &lowest_level_block_length_, sizeof(int));
    in.read((char *) &number_of_levels_, sizeof(int));
    in.read((char *) &rank_select_support_, sizeof(bool));

    for (int i = 0; i < number_of_levels_ - 1; ++i) {
        is_internal_.push_back(new sdsl::bit_vector());
        (*is_internal_[i]).load(in);
    }

    for (sdsl::bit_vector *bv : is_internal_) {
        bt_bv_rank_.push_back(new sdsl::rank_support_v<1>(bv));
    }

    for (int i = 0; i < number_of_levels_ - 1; ++i) {
        bt_offsets_.push_back(new sdsl::int_vector<>());
        (*bt_offsets_[i]).load(in);
    }

    leaf_string_ = new sdsl::int_vector<>();
    (*leaf_string_).load(in);

    alphabet_ = new sdsl::int_vector<>();
    (*alphabet_).load(in);

    int c = 0;
    for (int character : (*alphabet_)) {
        mapping_[character] = c++;
    }

    if (rank_select_support_) {
        for (int character : (*alphabet_)) {
            lowest_complete_level_prefix_ranks_[character] = new sdsl::int_vector<>();
            (*lowest_complete_level_prefix_ranks_[character]).load(in);
        }

        for (int character : (*alphabet_)) {
            for (int i = 0; i < number_of_levels_; ++i) {
                block_ranks_[character].push_back(new sdsl::int_vector<>());
                (*block_ranks_[character][i]).load(in);
            }
        }

        for (int character : (*alphabet_)) {
            for (int i = 0; i < number_of_levels_ - 1; ++i) {
                bt_second_ranks_[character].push_back(new sdsl::int_vector<>());
                (*bt_second_ranks_[character][i]).load(in);
            }
        }
    }
}

CBlockTree::~CBlockTree() {

    for (sdsl::bit_vector *bv : is_internal_) {
        delete bv;
    }

    for (sdsl::rank_support_v<1> *rank : bt_bv_rank_) {
        delete rank;
    }

    for (sdsl::int_vector<> *offsets : bt_offsets_) {
        delete offsets;
    }

    delete leaf_string_;
    delete alphabet_;

    if (rank_select_support_) {
        for (auto pair : lowest_complete_level_prefix_ranks_) {
            delete pair.second;
        }

        for (auto pair : block_ranks_) {
            for (sdsl::int_vector<> *ranks : pair.second) {
                delete ranks;
            }
        }

        for (auto pair : bt_second_ranks_) {
            for (sdsl::int_vector<> *ranks : pair.second) {
                delete ranks;
            }
        }
    }
}

int CBlockTree::access(int i) {

    int current_block  = i / lowest_level_block_length_;
    int current_length = lowest_level_block_length_;
    i -= current_block * lowest_level_block_length_;
    int level = 0;
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            current_length /= arity_;
            int child_number = i / current_length;
            i -= child_number * current_length;
            current_block = (*bt_bv_rank_[level])(current_block) *arity_ + child_number;
            ++level;
        } else { // Case BackBlock
            int encoded_offset = (*bt_offsets_[level])[current_block - (*bt_bv_rank_[level])(current_block + 1)];
            current_block      = encoded_offset / current_length;
            i += encoded_offset % current_length;
            if (i >= current_length) {
                ++current_block;
                i -= current_length;
            }
        }
    }
    return (*alphabet_)[(*leaf_string_)[i + current_block * current_length]];
}

int CBlockTree::rank(int c, int i) {

    auto &ranks        = block_ranks_[c];
    auto &second_ranks = bt_second_ranks_[c];

    int current_block  = i / lowest_level_block_length_;
    int current_length = lowest_level_block_length_;
    i                  = i - current_block * current_length;
    int level          = 0;

    int r = (*lowest_complete_level_prefix_ranks_[c])[current_block];
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            current_length /= arity_;
            int child_number = i / current_length;
            i -= child_number * current_length;

            int firstChild = (*bt_bv_rank_[level])(current_block) *arity_;
            for (int child = firstChild; child < firstChild + child_number; ++child) r += (*ranks[level + 1])[child];
            current_block = firstChild + child_number;
            ++level;
        } else { // Case BackBlock
            int index          = current_block - (*bt_bv_rank_[level])(current_block + 1);
            int encoded_offset = (*bt_offsets_[level])[index];
            current_block      = encoded_offset / current_length;
            i += encoded_offset % current_length;
            r += (*second_ranks[level])[index];
            if (i >= current_length) {
                ++current_block;
                i -= current_length;
            } else {
                r -= (*ranks[level])[current_block];
            }
        }
    }

    i += current_block * current_length;
    int d = mapping_[c];
    for (int j = current_block * current_length; j <= i; ++j) {
        if ((*leaf_string_)[j] == d)
            ++r;
    }

    return r;
}

int CBlockTree::select(int c, int k) {

    auto &ranks                    = block_ranks_[c];
    auto &second_ranks             = bt_second_ranks_[c];
    auto &first_level_prefix_ranks = lowest_complete_level_prefix_ranks_[c];

    int current_block = (k - 1) / lowest_level_block_length_;

    int end_block = (*first_level_prefix_ranks).size() - 1;
    while (current_block != end_block) {
        int m = current_block + (end_block - current_block) / 2;
        int f = (*first_level_prefix_ranks)[m];
        if (f < k) {
            if (end_block - current_block == 1) {
                if ((*first_level_prefix_ranks)[m + 1] < k) {
                    current_block = m + 1;
                }
                break;
            }
            current_block = m;
        } else {
            end_block = m - 1;
        }
    }

    int current_length = lowest_level_block_length_;
    int s              = current_block * current_length;
    k -= (*first_level_prefix_ranks)[current_block];
    int level = 0;
    while (level < number_of_levels_ - 1) {
        if ((*is_internal_[level])[current_block]) { // Case InternalBlock
            int firstChild          = (*bt_bv_rank_[level])(current_block) *arity_;
            int child               = firstChild;
            int r                   = (*ranks[level + 1])[child];
            int last_possible_child = (firstChild + arity_ - 1 > (*ranks[level + 1]).size() - 1)
                                          ? (*ranks[level + 1]).size() - 1
                                          : firstChild + arity_ - 1;
            while (child < last_possible_child && k > r) { // Border conditions?
                ++child;
                r += (*ranks[level + 1])[child];
            }
            k -= r - (*ranks[level + 1])[child];
            current_length /= arity_;
            s += (child - firstChild) * current_length;
            current_block = child;
            ++level;
        } else { // Case BackBlock
            int index          = current_block - (*bt_bv_rank_[level])(current_block + 1);
            int encoded_offset = (*bt_offsets_[level])[index];
            current_block      = encoded_offset / current_length;

            k -= (*second_ranks[level])[index];
            s -= encoded_offset % current_length;
            if (k > 0) {
                s += current_length;
                ++current_block;
            } else {
                k += (*ranks[level])[current_block];
            }
        }
    }

    int d = mapping_[c];
    for (int j = current_block * current_length;; ++j) {
        if ((*leaf_string_)[j] == d)
            --k;
        if (!k)
            return s + j - current_block * current_length;
    }

    return -1;
}

int CBlockTree::get_partial_size() {
    int fields = sizeof(int) * 3;

    int leaf_string_size = sdsl::size_in_bytes(*leaf_string_);

    int alphabet_size = sdsl::size_in_bytes(*alphabet_);
    int mapping_size  = sizeof(int) * 256;

    int bt_bv_size = sizeof(void *);
    for (sdsl::bit_vector *bv : is_internal_) {
        bt_bv_size += sdsl::size_in_bytes(*bv);
    }

    int bt_bv_rank_size = sizeof(void *);
    for (sdsl::rank_support_v<1> *bvr : bt_bv_rank_) {
        bt_bv_rank_size += sdsl::size_in_bytes(*bvr);
    }

    int bt_offsets_size = sizeof(void *);
    for (sdsl::int_vector<> *offsets : bt_offsets_) {
        bt_offsets_size += sdsl::size_in_bytes(*offsets);
    }

    return fields + mapping_size + alphabet_size + bt_bv_size + bt_bv_rank_size + bt_offsets_size + leaf_string_size;
}

int CBlockTree::size() {

    int bt_ranks_total_size = (block_ranks_.size() + 1) * sizeof(void *);
    for (auto pair : block_ranks_) {
        int size = 0;
        for (sdsl::int_vector<> *ranks : pair.second) {
            size += sdsl::size_in_bytes(*ranks);
        }
        bt_ranks_total_size += size;
    }

    int bt_prefix_ranks_first_level_size = 0;
    for (auto pair : lowest_complete_level_prefix_ranks_) {
        bt_prefix_ranks_first_level_size += sdsl::size_in_bytes(*(pair.second));
    }

    int bt_second_ranks_total_size = (bt_second_ranks_.size() + 1) * sizeof(void *);
    for (auto pair : bt_second_ranks_) {
        int size = 0;
        for (sdsl::int_vector<> *ranks : pair.second) {
            size += sdsl::size_in_bytes(*ranks);
        }
        bt_second_ranks_total_size += size;
    }

    int partial_total_size = get_partial_size();
    int rank_size          = bt_second_ranks_total_size + bt_ranks_total_size + bt_prefix_ranks_first_level_size;
    return rank_size + partial_total_size;
}

void CBlockTree::serialize(std::ostream &out) {

    out.write((char *) &arity_, sizeof(int));
    out.write((char *) &root_arity_, sizeof(int));
    out.write((char *) &lowest_level_block_length_, sizeof(int));
    out.write((char *) &number_of_levels_, sizeof(int));
    out.write((char *) &rank_select_support_, sizeof(bool));

    for (sdsl::bit_vector *bv : is_internal_) {
        (*bv).serialize(out);
    }

    for (sdsl::int_vector<> *offsets : bt_offsets_) {
        (*offsets).serialize(out);
    }

    (*leaf_string_).serialize(out);

    (*alphabet_).serialize(out);

    if (rank_select_support_) {
        for (int character : (*alphabet_)) {
            (*lowest_complete_level_prefix_ranks_[character]).serialize(out);
        }

        for (int character : (*alphabet_)) {
            for (sdsl::int_vector<> *ranks : block_ranks_[character]) {
                (*ranks).serialize(out);
            }
        }

        for (int character : (*alphabet_)) {
            for (sdsl::int_vector<> *second_ranks : bt_second_ranks_[character]) {
                (*second_ranks).serialize(out);
            }
        }
    }
}
