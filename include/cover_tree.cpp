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

