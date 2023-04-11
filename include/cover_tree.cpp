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

/****************************** k-Nearest Neighbours *************************************/


/****************************** Range Neighbours Search *************************************/


/****************************** Cover Trees Properties *************************************/





/****************************** Unit Tests for Cover Trees *************************************/

/****************************** Internal Constructors of Cover Trees *************************************/


/****************************** Public API for creation of Cover Trees *************************************/



/******************************************* Auxiliary Functions ***************************************************/



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



/******************************************* Functions to remove ***************************************************/

