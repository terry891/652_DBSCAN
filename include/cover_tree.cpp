/*
 * Copyright (c) 2017 Manzil Zaheer All rights reserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include "cover_tree.h"
#include "utils.h"

#include <numeric>

/******************************* Insert ***********************************************/

/******************************* Remove ***********************************************/


bool CoverTree::remove(const pointType &p)
{
    bool ret_val = false;
    // First find the point
    std::pair<CoverTree::Node*, double> result(root, root->dist(p));
    NearestNeighbour(root, result.second, p, result);

    if (result.second<=0.0)
    {   // point found
        CoverTree::Node* node_p = result.first;
        CoverTree::Node* parent_p = node_p->parent;
        if (node_p == root)
        {
            std::cout << "Sorry can not delete root efficiently!" << std::endl;
        }
        else
        {
            // 1. Remove p from parent's list of child
            unsigned num_children = parent_p->children.size();
            for (unsigned i = 0; i < num_children; ++i)
            {
                if (parent_p->children[i]==node_p)
                {
                    parent_p->children[i] =  parent_p->children.back();
                    parent_p->children.pop_back();
                    break;
                }

            }

            // 2. For each child q of p:
            for(CoverTree::Node* q : *node_p)
            {
                CoverTree::insert(root, q);
            }

            //3. delete
            delete node_p;

            ret_val = true;
        }
    }
    //calc_maxdist();
    return ret_val;
}


/****************************** Nearest Neighbour *************************************/
void CoverTree::NearestNeighbour(CoverTree::Node* current, double dist_current, const pointType &p, std::pair<CoverTree::Node*, double>& nn) const
{
    // If the current node is the nearest neighbour
    if (dist_current < nn.second)
    {
        nn.first = current;
        nn.second = dist_current;
    }

    // Sort the children
    unsigned num_children = current->children.size();
    std::vector<int> idx(num_children);
    std::iota(std::begin(idx), std::end(idx), 0);
    std::vector<double> dists(num_children);
    //dist_count[current->level].fetch_add(num_children, std::memory_order_relaxed);
    for (unsigned i = 0; i < num_children; ++i)
        dists[i] = current->children[i]->dist(p);
    auto comp_x = [&dists](int a, int b) { return dists[a] < dists[b]; };
    std::sort(std::begin(idx), std::end(idx), comp_x);

    for (const auto& child_idx : idx)
    {
        Node* child = current->children[child_idx];
        double dist_child = dists[child_idx];
        if (child->maxdistUB > current->covdist()/(base-1))
            std::cout << "I am crazy because max upper bound is bigger than 2**i " << child->maxdistUB << " " << current->covdist()/(base-1) << std::endl;
        if (nn.second > dist_child - child->maxdistUB)
            NearestNeighbour(child, dist_child, p, nn);
    }
}

// First the number of nearest neighbor
std::pair<CoverTree::Node*, double> CoverTree::NearestNeighbour(const pointType &p) const
{
    std::pair<CoverTree::Node*, double> result(root, root->dist(p));
    NearestNeighbour(root, result.second, p, result);
    return result;
}

/****************************** k-Nearest Neighbours *************************************/

void CoverTree::kNearestNeighbours(CoverTree::Node* current, double dist_current, const pointType& p, std::vector<std::pair<CoverTree::Node*, double>>& nnList) const
{
    // TODO(manzilz): An efficient implementation ?

    // If the current node is eligible to get into the list
    if(dist_current < nnList.back().second)
    {
        auto comp_x = [](std::pair<CoverTree::Node*, double> a, std::pair<CoverTree::Node*, double> b) { return a.second < b.second; };
        std::pair<CoverTree::Node*, double> temp(current, dist_current);
        nnList.insert(
            std::upper_bound( nnList.begin(), nnList.end(), temp, comp_x ),
            temp
        );
        nnList.pop_back();
    }

    // Sort the children
    unsigned num_children = current->children.size();
    std::vector<int> idx(num_children);
    std::iota(std::begin(idx), std::end(idx), 0);
    std::vector<double> dists(num_children);
    for (unsigned i = 0; i < num_children; ++i)
        dists[i] = current->children[i]->dist(p);
    auto comp_x = [&dists](int a, int b) { return dists[a] < dists[b]; };
    std::sort(std::begin(idx), std::end(idx), comp_x);

    for (const auto& child_idx : idx)
    {
        Node* child = current->children[child_idx];
        double dist_child = dists[child_idx];
        if ( nnList.back().second > dist_child - child->maxdistUB)
            kNearestNeighbours(child, dist_child, p, nnList);
    }
}

std::vector<std::pair<CoverTree::Node*, double>> CoverTree::kNearestNeighbours(const pointType &queryPt, unsigned numNbrs) const
{
    // Do the worst initialization
    std::pair<CoverTree::Node*, double> dummy(new CoverTree::Node(), std::numeric_limits<double>::max());
    // List of k-nearest points till now
    std::vector<std::pair<CoverTree::Node*, double>> nnList(numNbrs, dummy);

    // Call with root
    double dist_root = root->dist(queryPt);
    kNearestNeighbours(root, dist_root, queryPt, nnList);

    return nnList;
}

/****************************** Range Neighbours Search *************************************/

void CoverTree::rangeNeighbours(CoverTree::Node* current, double dist_current, const pointType &p, double range, std::vector<std::pair<CoverTree::Node*, double>>& nnList) const
{
    // If the current node is eligible to get into the list
    if (dist_current < range)
    {
        std::pair<CoverTree::Node*, double> temp(current, dist_current);
        nnList.push_back(temp);
    }

    // Sort the children
    unsigned num_children = current->children.size();
    std::vector<int> idx(num_children);
    std::iota(std::begin(idx), std::end(idx), 0);
    std::vector<double> dists(num_children);
    for (unsigned i = 0; i < num_children; ++i)
        dists[i] = current->children[i]->dist(p);
    auto comp_x = [&dists](int a, int b) { return dists[a] < dists[b]; };
    std::sort(std::begin(idx), std::end(idx), comp_x);

    for (const auto& child_idx : idx)
    {
        Node* child = current->children[child_idx];
        double dist_child = dists[child_idx];
        if (range > dist_child - child->maxdistUB)
            rangeNeighbours(child, dist_child, p, range, nnList);
    }
}

std::vector<std::pair<CoverTree::Node*, double>> CoverTree::rangeNeighbours(const pointType &queryPt, double range) const
{
    // List of nearest neighbors in the range
    std::vector<std::pair<CoverTree::Node*, double>> nnList;

    // Call with root
    double dist_root = root->dist(queryPt);
    rangeNeighbours(root, dist_root, queryPt, range, nnList);

    return nnList;
}

/****************************** Cover Trees Properties *************************************/

void CoverTree::generate_id(CoverTree::Node* current)
{
    // assign current node
    current->ID = N++;
#ifdef DEBUG
    std::cout << "Pre: " << current->ID << std::endl;
#endif

    // traverse children
    for (const auto& child : *current)
        generate_id(child);
}

// find true maxdist
void CoverTree::calc_maxdist()
{
    std::vector<CoverTree::Node*> travel;
    std::vector<CoverTree::Node*> active;

    CoverTree::Node* current = root;

    root->maxdistUB = 0.0;
    travel.push_back(root);
    while (!travel.empty())
    {
        current = travel.back();
        if (current->maxdistUB <= 0) {
            while (current->children.size()>0)
            {
                active.push_back(current);
                // push the children
                for (int i = current->children.size() - 1; i >= 0; --i)
                {
                    current->children[i]->maxdistUB = 0.0;
                    travel.push_back(current->children[i]);
                }
                current = current->children[0];
            }
        }
        else
            active.pop_back();

        // find distance with current node
        for (const auto& n : active)
            n->maxdistUB = std::max(n->maxdistUB, n->dist(current));

        // Pop
        travel.pop_back();
    }
}

/****************************** Serialization of Cover Trees *************************************/

// Pre-order traversal
char* CoverTree::preorder_pack(char* buff, CoverTree::Node* current) const
{
    // copy current node
    unsigned shift = current->_p.rows() * sizeof(pointType::Scalar);
    char* start = (char*)current->_p.data();
    char* end = start + shift;
    std::copy(start, end, buff);
    buff += shift;

    shift = sizeof(int);
    start = (char*)&(current->level);
    end = start + shift;
    std::copy(start, end, buff);
    buff += shift;

#ifdef DEBUG
    std::cout << "Pre: " << current->ID << std::endl;
#endif

    // traverse children
    for (const auto& child : *current)
        buff = preorder_pack(buff, child);

    return buff;
}

// Post-order traversal
char* CoverTree::postorder_pack(char* buff, CoverTree::Node* current) const
{
    // traverse children
    for (const auto& child : *current)
        buff = postorder_pack(buff, child);

    // save current node ID
#ifdef DEBUG
    std::cout << "Post: " << current->ID << std::endl;
#endif
    unsigned shift = sizeof(unsigned);
    char* start = (char*)&(current->ID);
    char* end = start + shift;
    std::copy(start, end, buff);
    buff += shift;

    return buff;
}

// reconstruct tree from Pre&Post traversals
void CoverTree::PrePost(CoverTree::Node*& current, char*& pre, char*& post)
{
    // The top element in pre-order list PRE is the root of T
    current = new CoverTree::Node;
    current->_p = pointType(D);
    for (unsigned i = 0; i < D; ++i)
    {
        current->_p[i] = *((pointType::Scalar *)pre);
        pre += sizeof(pointType::Scalar);
    }
    current->level = *((int *)pre);
    current->ID = N++;
    current->maxdistUB = 0;
    pre += sizeof(int);

    // Construct sub-trees until the root is found in the post-order list
    while (*((unsigned*)post) != current->ID)
    {
        CoverTree::Node* temp = NULL;
        PrePost(temp, pre, post);
        current->children.push_back(temp);
    }

    // All sub-trees of T are constructed
    post += sizeof(unsigned); // Delete top element of POST
}

unsigned CoverTree::msg_size() const
{
    return 2 * sizeof(unsigned)
        + sizeof(pointType::Scalar)*D*N
        + sizeof(int)*N
        + sizeof(unsigned)*N;
}

// Serialize to a buffer
char* CoverTree::serialize() const
{
    //// check if valid id present
    //if (!id_valid)
    //{
    // N = 0;
    // generate_id(root);
    // id_valid = true;
    //}
    //count_points();

    // Covert following to char* buff with following order
    // N | D | (points, levels) | List
    char* buff = new char[msg_size()];

    char* pos = buff;

    // insert N
    unsigned shift = sizeof(unsigned);
    char* start = (char*)&(N);
    char* end = start + shift;
    std::copy(start, end, pos);
    pos += shift;

    // insert D
    shift = sizeof(unsigned);
    start = (char*)&(D);
    end = start + shift;
    std::copy(start, end, pos);
    pos += shift;

    // insert points and level
    pos = preorder_pack(pos, root);
    pos = postorder_pack(pos, root);

    return buff;
}

// Deserialize from a buffer
void CoverTree::deserialize(char* buff)
{
    /** Convert char* buff into following buff = N | D | (points, levels) | List **/
    char* save = buff;

    // extract N and D
    N = *((unsigned *)buff);
    buff += sizeof(unsigned);
    D = *((unsigned *)buff);
    buff += sizeof(unsigned);

    // pointer to post-order list
    char* post = buff + sizeof(pointType::Scalar)*D*N
        + sizeof(int)*N;

    // reconstruction
    N = 0;
    PrePost(root, buff, post);

    delete[] save;
}

/****************************** Unit Tests for Cover Trees *************************************/
bool CoverTree::check_covering() const
{
    bool result = true;
    std::stack<CoverTree::Node*> travel;
    CoverTree::Node* curNode;

    // Initialize with root
    travel.push(root);

    // Pop, check and then push the children
    while (!travel.empty())
    {
        // Pop
        curNode = travel.top();
        travel.pop();

        // Check covering for the current -> children pair
        for (const auto& child : *curNode)
        {
            travel.push(child);
            if( curNode->dist(child) > curNode->covdist() )
                result = false;
            //std::cout << *curNode << " -> " << *child << " @ " << curNode->dist(child) << " | " << curNode->covdist() << std::endl;
        }
    }

    return result;
}

/****************************** Internal Constructors of Cover Trees *************************************/

// constructor: NULL tree
CoverTree::CoverTree(const int truncate /*=-1*/ )
    : root(NULL)
    , min_scale(1000)
    , max_scale(0)
    , truncate_level(truncate)
    , id_valid(false)
    , N(0)
    , D(0)
{
}

// constructor: needs at least 1 point to make a valid cover-tree
CoverTree::CoverTree(const pointType& p, int truncate /*=-1*/)
    : min_scale(1000)
    , max_scale(0)
    , truncate_level(truncate)
    , id_valid(false)
    , N(1)
    , D(p.rows())
{
    root = new CoverTree::Node;
    root->_p = p;
    root->level = 0;
    root->maxdistUB = 0;
}


/****************************** Public API for creation of Cover Trees *************************************/



/******************************************* Auxiliary Functions ***************************************************/

// get root level == max_level
int CoverTree::get_level()
{
    return root->level;
}

void CoverTree::print_levels()
{
    std::stack<CoverTree::Node*> travel;
    CoverTree::Node* curNode;

    // Initialize with root
    travel.push(root);

    // Pop, print and then push the children
    while (!travel.empty())
    {
        // Pop
        curNode = travel.top();
        travel.pop();

        // Count the level
        level_count[curNode->level]++;

        // Now push the children
        for (const auto& child : *curNode)
            travel.push(child);
    }

    for(auto const& qc : level_count)
    {
        std::cout << "Number of nodes at level " << qc.first << " = " << qc.second << std::endl;
        //dist_count[qc.first].store(0);
    }
}

// Pretty print
/*
std::ostream& operator<<(std::ostream& os, const CoverTree& ct)
{
    std::stack<CoverTree::Node*> travel;
    CoverTree::Node* curNode;

    // Initialize with root
    travel.push(ct.root);

    // Qualitatively keep track of number of prints
    int numPrints = 0;
    // Pop, print and then push the children
    while (!travel.empty())
    {
        if (numPrints > 5000)
            throw std::runtime_error("Printing stopped prematurely, something wrong!");
        numPrints++;

        // Pop
        curNode = travel.top();
        travel.pop();

        // Print the current -> children pair
        for (const auto& child : *curNode)
            os << *curNode << " -> " << *child << std::endl;

        // Now push the children
        for (int i = curNode->children.size() - 1; i >= 0; --i)
            travel.push(curNode->children[i]);
    }

    return os;
}
*/

std::vector<pointType> CoverTree::get_points()
{
    std::vector<pointType> points;

    std::stack<CoverTree::Node*> travel;
    CoverTree::Node* current;

    // Initialize with root
    travel.push(root);

    N = 0;
    // Pop, print and then push the children
    while (!travel.empty())
    {
        // Pop
        current = travel.top();
        travel.pop();

        // Add to dataset
        points.push_back(current->_p);
        current->ID = N++;

        // Now push the children
        for (const auto& child : *current)
            travel.push(child);
    }

    return points;
}


/******************************************* Functions to remove ***************************************************/

unsigned CoverTree::count_points()
{
    std::stack<CoverTree::Node*> travel;
    CoverTree::Node* current;

    // Initialize with root
    travel.push(root);

    unsigned result = 0;
    // Pop, print and then push the children
    while (!travel.empty())
    {
        // Pop
        current = travel.top();
        travel.pop();

        // Add to dataset
        ++result;

        // Now push the children
        for (const auto& child : *current)
            travel.push(child);
    }
    return result;
}

