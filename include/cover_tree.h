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

# ifndef _COVER_TREE_H
# define _COVER_TREE_H

//#define DEBUG

#include <atomic>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <vector>
#include <shared_mutex>

#include "dataset.h"
#include "utils.h"

#ifdef __clang__
#define SHARED_MUTEX_TYPE shared_mutex
#else
#define SHARED_MUTEX_TYPE shared_timed_mutex
#endif

#include <Eigen/Core>
typedef Eigen::VectorXd pointType;
//typedef pointType::Scalar dtype;

template <typename T>
class CoverTree
{
/************************* Internal Functions ***********************************************/
public:
    /*** Base to use for the calculations ***/
    static constexpr double base = 1.3;
    static constexpr int SIZE=2048;
    //std::array<double,2048> compute_pow_table(){
    //  std::array<double,2048> powdict{};
    //  powdict[0] = 2.098975e-117;
    //  for (int i = 1; i<2048; ++i)
    //      powdict[i] = 1.3* powdict[i-1];
    //  return powdict;
    //}
    static double powdict[SIZE];

public:
    /*** structure for each node ***/
    struct Node
    {
        pointType _p;                       // point associated with the node
        std::vector<Node*> children;        // list of children
        int level;                          // current level of the node
        double maxdistUB;                   // upper bound of distance to any of descendants
        unsigned ID;                        // unique ID of current node
        Node* parent;                       // parent of current node

        mutable std::SHARED_MUTEX_TYPE mut;// lock for current node

        /*** Node modifiers ***/
        double covdist()                    // covering distance of subtree at current node
        {
            return powdict[level + 1024];
        }
        double sepdist()                    // separating distance between nodes at current level
        {
            return powdict[level + 1023];
        }
        double dist(const pointType& pp) const  // L2 distance between current node and point pp
        {
            return (_p - pp).norm();
        }
        double dist(Node* n) const              // L2 distance between current node and node n
        {
            return (_p - n->_p).norm();
        }
        Node* setChild(const pointType& pIns, int new_id=-1)   // insert a new child of current node with point pIns
        {
            Node* temp = new Node;
            temp->_p = pIns;
            temp->level = level - 1;
            temp->maxdistUB = 0; // powdict[level + 1024];
            temp->ID = new_id;
            temp->parent = this;
            children.push_back(temp);
            return temp;
        }
        Node* setChild(Node* pIns)          // insert the subtree pIns as child of current node
        {
            if( pIns->level != level - 1)
            {
                Node* current = pIns;
                std::stack<Node*> travel;
                current->level = level-1;
                //current->maxdistUB = powdict[level + 1024];
                travel.push(current);
                while (!travel.empty())
                {
                    current = travel.top();
                    travel.pop();

                    for (const auto& child : *current)
                    {
                        child->level = current->level-1;
                        //child->maxdistUB = powdict[child->level + 1025];
                        travel.push(child);
                    }

                }
            }
            pIns->parent = this;
            children.push_back(pIns);
            return pIns;
        }

        /*** erase child ***/
        void erase(size_t pos)
        {
            children[pos] = children.back();
            children.pop_back();
        }

//        void erase(std::vector<CoverTree<T>::Node*>::iterator pos)
//        {
//            *pos = children.back();
//            children.pop_back();
//        }
//
//        /*** Iterator access ***/
//        inline std::vector<CoverTree<T>::Node*>::iterator begin()
//        {
//            return children.begin();
//        }
//        inline std::vector<CoverTree<T>::Node*>::iterator end()
//        {
//            return children.end();
//        }
//        inline std::vector<CoverTree<T>::Node*>::const_iterator begin() const
//        {
//            return children.begin();
//        }
//        inline std::vector<CoverTree<T>::Node*>::const_iterator end() const
//        {
//            return children.end();
//        }
    };
    // mutable std::map<int,std::atomic<unsigned>> dist_count;
    std::map<int,unsigned> level_count;

protected:
    Node* root;                         // Root of the tree
    std::atomic<int> min_scale;         // Minimum scale
    std::atomic<int> max_scale;         // Minimum scale
    //int min_scale;                    // Minimum scale
    //int max_scale;                    // Minimum scale
    int truncate_level;                 // Relative level below which the tree is truncated
    bool id_valid;
    std::vector<T> mean_data;

    std::atomic<size_t> N;            // Number of points in the cover tree
    //unsigned N;                       // Number of points in the cover tree
    unsigned D;                         // Dimension of the points

    std::SHARED_MUTEX_TYPE global_mut;  // lock for changing the root

    /*** Insert point or node at current node ***/
    bool insert(Node* current, const pointType& p){
      bool result = false;
      if (truncate_level > 0 && current->level < max_scale-truncate_level)
          return false;
  
      // acquire read lock
      current->mut.lock_shared();
  
      // Sort the children
      unsigned num_children = current->children.size();
      std::vector<int> idx(num_children);
      std::iota(std::begin(idx), std::end(idx), 0);
      std::vector<double> dists(num_children);
      for (unsigned i = 0; i < num_children; ++i)
          dists[i] = current->children[i]->dist(p);
      auto comp_x = [&dists](int a, int b) { return dists[a] < dists[b]; };
      std::sort(std::begin(idx), std::end(idx), comp_x);
  
      bool flag = true;
      for (const auto& child_idx : idx)
      {
          Node* child = current->children[child_idx];
          double dist_child = dists[child_idx];
          if (dist_child <= 0.0)
          {
              // release read lock then enter child
              current->mut.unlock_shared();
              flag = false;
              std::cout << "Duplicate entry!!!" << std::endl;
              break;
          }
          else if (dist_child <= child->covdist())
          {
              // release read lock then enter child
              if (child->maxdistUB < dist_child)
                  child->maxdistUB = dist_child;
              current->mut.unlock_shared();
              result = insert(child, p);
              flag = false;
              break;
          }
      }
  
      if (flag)
      {
          // release read lock then acquire write lock
          current->mut.unlock_shared();
          current->mut.lock();
          // check if insert is still valid, i.e. no other point was inserted else restart
          if (num_children==current->children.size())
          {
              int new_id = ++N;
              current->setChild(p, new_id);
              result = true;
              current->mut.unlock();
  
              int local_min = min_scale.load();
              while( local_min > current->level - 1){
                  min_scale.compare_exchange_weak(local_min, current->level - 1, std::memory_order_relaxed, std::memory_order_relaxed);
                  local_min = min_scale.load();
              }
          }
          else
          {
              current->mut.unlock();
              result = insert(current, p);
          }
      }
      return result;

    }
    bool insert(Node* current, Node* p){
      bool result = false;
      std::cout << "Node insert called!";
      if (truncate_level > 0 && current->level < max_scale-truncate_level)
          return false;

      // acquire read lock
      current->mut.lock_shared();

      // Sort the children
      unsigned num_children = current->children.size();
      std::vector<int> idx(num_children);
      std::iota(std::begin(idx), std::end(idx), 0);
      std::vector<double> dists(num_children);
      for (unsigned i = 0; i < num_children; ++i)
          dists[i] = current->children[i]->dist(p);
      auto comp_x = [&dists](int a, int b) { return dists[a] < dists[b]; };
      std::sort(std::begin(idx), std::end(idx), comp_x);

      bool flag = true;
      for (const auto& child_idx : idx)
      {
          Node* child = current->children[child_idx];
          double dist_child = dists[child_idx];
          if (dist_child <= 0.0)
          {
              // release read lock then enter child
              current->mut.unlock_shared();
              flag = false;
              break;
          }
          else if (dist_child <= child->covdist())
          {
              // release read lock then enter child
              current->mut.unlock_shared();
              result = insert(child, p);
              flag = false;
              break;
          }
      }

      if (flag)
      {
          // release read lock then acquire write lock
          current->mut.unlock_shared();
          current->mut.lock();
          // check if insert is still valid, i.e. no other point was inserted else restart
          if (num_children==current->children.size())
          {
              ++N;
              current->setChild(p);
              result = true;
              current->mut.unlock();

              int local_min = min_scale.load();
              while( local_min > current->level - 1){
                  min_scale.compare_exchange_weak(local_min, current->level - 1, std::memory_order_relaxed, std::memory_order_relaxed);
                  local_min = min_scale.load();
              }
          }
          else
          {
              current->mut.unlock();
              result = insert(current, p);
          }
      }
      return result;

    }

    /*** Nearest Neighbour search ***/
    void NearestNeighbour(Node* current, double dist_current, const pointType &p, std::pair<CoverTree::Node*, double>& nn) const;

    /*** k-Nearest Neighbour search ***/
    void kNearestNeighbours(Node* current, double dist_current, const pointType& p, std::vector<std::pair<CoverTree::Node*, double>>& nnList) const;

    /*** Range search ***/
    void rangeNeighbours(Node* current, double dist_current, const pointType &p, double range, std::vector<std::pair<CoverTree::Node*, double>>& nnList) const;

    /*** Serialize/Desrialize helper function ***/
    char* preorder_pack(char* buff, Node* current) const;       // Pre-order traversal
    char* postorder_pack(char* buff, Node* current) const;      // Post-order traversal
    void PrePost(Node*& current, char*& pre, char*& post);

    /*** debug functions ***/
    unsigned msg_size() const;
    void calc_maxdist();                            //find true maxdist
    void generate_id(Node* current);                //Generate IDs for each node from root as 0


/************************* Public API ***********************************************/
public:
    CoverTree(Dataset& m_data){
      //CoverTree<T>::powdict  = CoverTree<T>::compute_pow_table();
      for(int i=0;i<SIZE;i++){
        powdict[i] = std::pow(1.3, i-1024);
      }
      std::vector<pointType> pList = from_dataset(m_data);
      int truncate = -1;
      int end = pList.size();
      int begin = 0;

      //1. Compute the mean of entire data
      pointType mx = utils::ParallelAddList(pList).get_result()/pList.size();

      //2. Compute distance of every point from the mean || Variance
      pointType dists = utils::ParallelDistanceComputeList(pList, mx).get_result();

      //3. argort the distance to find approximate mediod
      std::vector<int> idx(end-begin);
      std::iota(std::begin(idx), std::end(idx), 0);
      auto comp_x = [&dists](int a, int b) { return dists[a] > dists[b]; };
      std::sort(std::begin(idx), std::end(idx), comp_x);
      std::cout<<"Max distance: " << dists[idx[0]] << std::endl;

      //4. Compute distance of every point from the mediod
      mx = pList[idx[0]];
      dists = utils::ParallelDistanceComputeList(pList, mx).get_result();

      int scale_val = std::ceil(std::log(dists.maxCoeff())/std::log(base));
      std::cout<<"Scale chosen: " << scale_val << std::endl;
      pointType temp = pList[idx[0]];
      min_scale = scale_val; //-1000;
      max_scale = scale_val; //-1000;
      truncate_level = truncate;
      N = 1;
      D = temp.rows();

/*
      root = new CoverTree::Node;
      root->_p = temp;
      root->level = scale_val; //-1000;
      root->maxdistUB = powdict[scale_val+1024];
      int run_till = 50000 < end ? 50000 : end;
      for (int i = 1; i < run_till; ++i){
          //utils::progressbar(i, run_till);
          if(!insert(pList[idx[i]]))
              std::cout << "Insert failed!!!" << std::endl;
      }
      utils::progressbar(run_till, run_till);
      std::cout<<std::endl;

      std::cout << pList[0].rows() << ", " << pList.size() << std::endl;

      utils::parallel_for_progressbar(50000,end,[&](int i)->void{
      //for (int i = 50000; i < end; ++i){
          //utils::progressbar(i, end-50000);
          if(!insert(pList[idx[i]]))
              std::cout << "Insert failed!!!" << std::endl;
      });
*/
    }

    explicit CoverTree(int truncate = -1);
    // cover tree with one point as root
    CoverTree(const pointType& p, int truncate = -1);

    /*** Destructor ***/
    /*** Destructor: deallocating all memories by a post order traversal ***/
    ~CoverTree(){}

    static std::vector<pointType> from_dataset(Dataset& m_data){
      const size_t items = m_data.m_chunk[0];
      const size_t dimensions = m_data.m_chunk[1];
    
      std::vector<pointType> pList(items);
      // pragma omp optimizations for parallel constructions
      for (size_t i=0;i<items;i++){
        const T* point = static_cast<T*>(m_data.m_p) + i * dimensions;
    
        pointType newPoint = pointType(dimensions);
        for(size_t d=0;d<dimensions;d++){
          newPoint[d] = point[d];
        }
        pList.push_back(newPoint);
      }
    
      return pList;
    }

    /*** Insert point p into the cover tree ***/
    bool insert(const pointType& p){
      bool result = false;
      id_valid = false;
      global_mut.lock_shared();
      if (root->dist(p) > root->covdist())
      {
          global_mut.unlock_shared();
          std::cout<<"Entered case 1: " << root->dist(p) << " " << root->covdist() << " " << root->level <<std::endl;
          std::cout<<"Requesting global lock!" <<std::endl;
          global_mut.lock();
          while (root->dist(p) > base * root->covdist()/(base-1))
          {
              CoverTree::Node* current = root;
              CoverTree::Node* parent = NULL;
              while (current->children.size()>0)
              {
                  parent = current;
                  current = current->children.back();
              }
              if (parent != NULL)
              {
                  parent->children.pop_back();
                  current->level = root->level + 1;
                  //current->maxdistUB = powdict[current->level + 1025];
                  current->children.push_back(root);
                  root = current;
              }
              else
              {
                  root->level += 1;
                  //root->maxdistUB = powdict[root->level + 1025];
              }
          }
          ++N;
          CoverTree::Node* temp = new CoverTree::Node;
          temp->_p = p;
          temp->level = root->level + 1;
          temp->parent = NULL;
          //temp->maxdistUB = powdict[temp->level+1025];
          temp->children.push_back(root);
          root->parent = temp;
          root = temp;
          max_scale = root->level;
          result = true;
          //std::cout << "Upward: " << minScale << " " << maxScale << std::endl;
          global_mut.unlock();
          global_mut.lock_shared();
      }
      else
      {
          //root->tempDist = root->dist(p);
          result = insert(root, p);
      }
      global_mut.unlock_shared();
      return result;

    }

    /*** Remove point p into the cover tree ***/
    bool remove(const pointType& p);

    /*** Nearest Neighbour search ***/
    std::pair<CoverTree::Node*, double> NearestNeighbour(const pointType &p) const;

    /*** k-Nearest Neighbour search ***/
    std::vector<std::pair<CoverTree::Node*, double>> kNearestNeighbours(const pointType &p, unsigned k = 10) const;

    /*** Range search ***/
    std::vector<std::pair<CoverTree::Node*, double>> rangeNeighbours(const pointType &queryPt, double range = 1.0) const;

    /*** Serialize/Desrialize: useful for MPI ***/
    char* serialize() const;                                    // Serialize to a buffer
    void deserialize(char* buff);                               // Deserialize from a buffer

    /*** Unit Tests ***/
    bool check_covering() const;

    /*** Return the level of root in the cover tree (== max_level) ***/
    int get_level();
    void print_levels();

    /*** Return all points in the tree ***/
    std::vector<pointType> get_points();

    /*** Count the points in the tree ***/
    unsigned count_points();

    /*** Pretty print ***/
    //friend std::ostream& operator<<(std::ostream& os, const CoverTree& ct);
};

#endif //_COVER_TREE_H

