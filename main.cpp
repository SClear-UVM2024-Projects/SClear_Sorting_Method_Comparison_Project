/**
 * Samuel Clear
 * 9/10/2021
 */
#include <iostream>
#include "LoanRecipient.h"
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>

// These two include statements will be used to import a clock function into the program. This
// imported file was references from the following link:
// https://stackoverflow.com/questions/3220477/how-to-use-clock-in-c
#include <cstdio>
#include <ctime>

using std::vector;
using std::ifstream;
using ::clock_t;

/**
 * Main Function
 */

// Heap Sort

template<typename Comparable>
void printVec(const vector<Comparable> &v) {
    for (int i = 0; i < v.size(); ++i) {
        if (i != 0) {
            cout << ", ";
        }
        cout << v[i]; // I have decided to remove the single read that prints from each vector,
        // since the vectors will not be printed in the final version of the main function file.
    }
    cout << endl;
}

// Helper function for heapSort
inline int leftChild(int i) {
    return 2 * i + 1;
}

// Helper function for heapSort
// i is the index of the value being percolated down
// n is the number of items in the heap
template <typename Comparable>
void percolateDown(vector<Comparable> &items, int i, int n, int child, Comparable tmp, int& reads, int& writes) {
    for(tmp = items[i]; leftChild(i) < n; i = child) { // 1 read, 1 write
        child = leftChild(i);
        // choose the child with the larger value
        if (child != n - 1 && items[child] < items[child + 1]) { // 2 reads
            ++child;
        }
        reads += 3;
        writes += 1;
        // if the parent is less than the child, swap them
        if (tmp < items[child]) { // 2 reads
            items[i] = items[child]; // One read, one write
            reads += 3;
            writes += 1;
        } else { // 2 reads
            // parent is >= both children. nothing more to do.
            reads += 2;
            break;
        }
    }
    items[i] = tmp; // 1 write and 1 read
    ++reads;
    ++writes;
}

template <typename Comparable>
void heapSort(vector<Comparable> items, int& reads, int& writes) {
    int i, j, child;
    Comparable temp, tmp;
    // build the heap (with max value at root)
    for (i = items.size() / 2 - 1; i >= 0; --i) {
        percolateDown(items, i, items.size(), child, tmp, reads, writes); // Reference Percolate Down for Reads and Writes
    }
    // printVec(items); // Reference PrintVec for Reads and Writes
    // keep deleting the max
    for (j = items.size() - 1; j > 0; --j) {
        // swap the maximum out
        temp = items[0];
        items[0] = items[j];
        items[j] = temp; // 3 reads and 3 writes
        reads += 3;
        writes += 3;
        // make it into a heap again
        percolateDown(items, 0, j, child, tmp, reads, writes);
        //printVec(items);
        ++reads;
    }
}

// Bubble Sort

template<typename Comparable>
void bubbleSort(vector<Comparable> vec, int& reads, int& writes) {
    bool haveSwapped = true;
    int maxIndex = vec.size(), i;
    Comparable temp;
    while (haveSwapped) {
        haveSwapped = false;
        for (i = 0; i + 1 < maxIndex; ++i) {
            // Compare items at indices i and i+1 and swap if necessary
            if (vec[i] > vec[i+1]) { // 2 reads
                temp = vec[i];
                vec[i] = vec[i+1];
                vec[i+1] = temp; // 3 read and 3 writes
                // Update haveSwapped
                haveSwapped = true;
                reads += 3;
                writes += 3;
            }
            reads += 2; // Items are compared whether or not the statement is true
        }
        // printVec(); // Reference PrintVec for Reads and Writes
        // Update maxIndex
        --maxIndex;
    }
}

// Selection Sort

template<typename Comparable>
void selectionSort(vector<Comparable> vec, int& reads, int& writes) {
    int swapIndex, i, minIndex;
    Comparable temp;
    for (swapIndex = 0; swapIndex < vec.size() - 1; ++swapIndex) {
        // Loop through vector starting at swapIndex and keep track of min
        minIndex = swapIndex;
        for (i = swapIndex + 1; i < vec.size(); ++i) {
            if (vec[i] < vec[minIndex]) { // 2 reads
                minIndex = i;
            }
            reads += 2;
        }
        //  Observation:: This sorting method has an extremely low number of reading and writing methods because
        //  the indexes are being recorded, rather than the actual objects themselves.
        // Swap min value into swapIndex spot in vector
        temp = vec[swapIndex];
        vec[swapIndex] = vec[minIndex];
        vec[minIndex] = temp; // 3 reads and 3 writes
        reads += 3;
        writes += 3;
    }
}

// Insertion Sort
template<typename Comparable>
void insertionSort(vector<Comparable> vec, int& reads, int& writes) {
    int unsortedStartIndex, insertIndex;
    Comparable toBeInserted;
    for (unsortedStartIndex = 1; unsortedStartIndex < vec.size(); ++unsortedStartIndex) {
        toBeInserted = vec[unsortedStartIndex]; // 1 read and 1 write
        // Loop to shift over the larger elements
        insertIndex = unsortedStartIndex - 1;
        //  This loop must make 2 reads to run through another iteration of the loop. Whether the
        //  loop iterates again or not, the two reads must be made even if the loop does not iterate.
        //  This is why the reads accumulator variable adds two after the while loop ends.
        while (insertIndex >= 0 && vec[insertIndex] > toBeInserted) { // 2 reads
            vec[insertIndex + 1] = vec[insertIndex]; // 1 read and 1 write
            --insertIndex;
            reads += 3; // +2 reads from each iteration of this loop
            ++writes;
        }
        reads += 2; // 2 reads after the loop is told to conclude
        // Put toBeInserted back into vec
        vec[insertIndex + 1] = toBeInserted; // 1 read and 1 write
        reads += 2;
        writes += 2;
    }
}

// Merge Sort

template<typename Comparable>
void mergeSortRec(vector<Comparable> &vec, int startIndex, int endIndex, int& reads, int& writes) {
    // Recursive base case
    if (startIndex == endIndex) {
        // We have one item. There is nothing to split or sort.
        return;
    }

    // Recursive calls
    int centerIndex = (startIndex + endIndex) / 2;
    // Each of the two calls change the total number of reads and writes of the
    // parameterized variable passed the function.
    mergeSortRec(vec, startIndex, centerIndex, reads, writes);
    mergeSortRec(vec, centerIndex + 1, endIndex, reads, writes);

    // Merge
    vector<Comparable> temp;
    int leftIndex = startIndex;
    int rightIndex = centerIndex + 1;
    while (leftIndex <= centerIndex && rightIndex <= endIndex) {
        if (vec[leftIndex] <= vec[rightIndex]) { // 2 reads
            temp.push_back(vec[leftIndex]); // 1 write a 1 read
            ++leftIndex;
            reads += 3;
            ++writes;
        } else { // 2 reads
            temp.push_back(vec[rightIndex]); // 1 write and 1 read
            ++rightIndex;
            ++writes;
            reads += 3;
        }
    }
    // At this point, one of the halves has been completely copied into temp but the other hasn't
    // We need to finish copying everything into temp, so we make loops for each half
    while (leftIndex <= centerIndex) {
        temp.push_back(vec[leftIndex]); // 1 write and 1 read
        ++leftIndex;
        ++writes;
        ++reads;
    }
    while (rightIndex <= endIndex) {
        temp.push_back(vec[rightIndex]); // 1 write and 1 read
        ++rightIndex;
        ++writes;
        ++reads;
    }
    // At this point, all of the items from startIndex to endIndex have been copied into temp
    // Copy everything from temp back into vec
    for (int i = 0; i < temp.size(); ++i) {
        vec[i + startIndex] = temp[i]; // 1 write and 1 read
        ++writes;
        ++reads;
    }
}

template<typename Comparable>
void mergeSort(vector<Comparable> vec, int& reads, int& writes) {
    mergeSortRec(vec, 0, vec.size() - 1, reads, writes);
}

// Two-Sort, with one Vector of Integers sorted using Bubble Sort and one vector of
// LoanRecipient getID unique integer values sorted using Selection Sort

// Each of the versions of Two-Sort sort the lists by the variable YoE in each of the
// LoanRecipient objects which contains duplicates. This variable is used to demonstrate
// the stability of the second sorting function of the Two-Sort method.
void bubbleSortTwoSort(vector<LoanRecipient> &vec, int& reads, int& writes) {
    bool haveSwapped = true;
    int maxIndex = vec.size(), i;
    LoanRecipient temp;
    while (haveSwapped) {
        haveSwapped = false;
        for (i = 0; i + 1 < maxIndex; ++i) {
            // Compare items at indices i and i+1 and swap if necessary
            if (vec[i].getYearsOfExp() > vec[i + 1].getYearsOfExp()) { // 2 reads
                temp = vec[i];
                vec[i] = vec[i + 1];
                vec[i + 1] = temp; // 3 read and 3 writes
                // Update haveSwapped
                haveSwapped = true;
                reads += 3;
                writes += 3;
            }
            reads += 2;
        }
        // Reference PrintVec for Reads and Writes
        // Update maxIndex
        --maxIndex;
    }
}

// This is the sorting method used for the first sort of the Two-Sort method for the first version of Two-Sort
// Two versions of this selection sort method are created to allow for the different
// second sorting functions to be called.
template<typename Comparable>
void selectionSortV1(vector<Comparable> vec, int& reads, int& writes) {
    int swapIndex, i, minIndex;
    Comparable temp;
    for (swapIndex = 0; swapIndex < vec.size() - 1; ++swapIndex) {
        // Loop through vector starting at swapIndex and keep track of min
        minIndex = swapIndex;
        for (i = swapIndex + 1; i < vec.size(); ++i) {
            if (vec[i] < vec[minIndex]) { // 2 reads
                minIndex = i;
            }
            reads += 2;
        }
        //  Observation:: This sorting method has an extremely low number of reading and writing methods because
        //  the indexes are being recorded, rather than the actual objects themselves.
        // Swap min value into swapIndex spot in vector
        temp = vec[swapIndex];
        vec[swapIndex] = vec[minIndex];
        vec[minIndex] = temp; // 3 reads and 3 writes
        reads += 3;
        writes += 3;
    }
    bubbleSortTwoSort(vec, reads, writes);
}

// This is the sorting method used for the first sort of the Two-Sort method for the second version of Two-Sort
template<typename Comparable>
void selectionSortV2(vector<Comparable> vec, int& reads, int& writes) {
    int swapIndex, i, minIndex;
    Comparable temp;
    for (swapIndex = 0; swapIndex < vec.size() - 1; ++swapIndex) {
        // Loop through vector starting at swapIndex and keep track of min
        minIndex = swapIndex;
        for (i = swapIndex + 1; i < vec.size(); ++i) {
            if (vec[i] < vec[minIndex]) { // 2 reads
                minIndex = i;
            }
            reads += 2;
        }
        //  Observation:: This sorting method has an extremely low number of reading and writing methods because
        //  the indexes are being recorded, rather than the actual objects themselves.
        // Swap min value into swapIndex spot in vector
        temp = vec[swapIndex];
        vec[swapIndex] = vec[minIndex];
        vec[minIndex] = temp; // 3 reads and 3 writes
        reads += 3;
        writes += 3;
    }
    mergeSortTwoSort(vec, reads, writes);
}

// This version of Merge Sort that is used in the custom sorting method, TwoSort.
void mergeSortTwoSortRec(vector<LoanRecipient> &vec, int startIndex, int endIndex, int& reads, int& writes) {
    // Recursive base case
    if (startIndex == endIndex) {
        // We have one item. There is nothing to split or sort.
        return;
    }

    // Recursive calls
    int centerIndex = (startIndex + endIndex) / 2;
    // Each of the two calls change the total number of reads and writes of the
    // parameterized variable passed the function.
    mergeSortTwoSortRec(vec, startIndex, centerIndex, reads, writes);
    mergeSortTwoSortRec(vec, centerIndex + 1, endIndex, reads, writes);

    // Merge
    vector<LoanRecipient> temp;
    int leftIndex = startIndex;
    int rightIndex = centerIndex + 1;
    while (leftIndex <= centerIndex && rightIndex <= endIndex) {
        // Each of the LoanRecipient items are sorted by their YoE value, which maintains
        // the same order for duplicate items that were sorted in the previous sorting algorithm
        if (vec[leftIndex].getYearsOfExp() <= vec[rightIndex].getYearsOfExp()) { // 2 reads
            temp.push_back(vec[leftIndex]); // 1 write a 1 read
            ++leftIndex;
            reads += 3;
            ++writes;
        } else { // 2 reads
            temp.push_back(vec[rightIndex]); // 1 write and 1 read
            ++rightIndex;
            ++writes;
            reads += 3;
        }
    }
    // At this point, one of the halves has been completely copied into temp but the other hasn't
    // We need to finish copying everything into temp, so we make loops for each half
    while (leftIndex <= centerIndex) {
        temp.push_back(vec[leftIndex]); // 1 write and 1 read
        ++leftIndex;
        ++writes;
        ++reads;
    }
    while (rightIndex <= endIndex) {
        temp.push_back(vec[rightIndex]); // 1 write and 1 read
        ++rightIndex;
        ++writes;
        ++reads;
    }
    // At this point, all the items from startIndex to endIndex have been copied into temp
    // Copy everything from temp back into vec
    for (int i = 0; i < temp.size(); ++i) {
        vec[i + startIndex] = temp[i]; // 1 write and 1 read
        ++writes;
        ++reads;
    }
    // printVec(vec);
}

void mergeSortTwoSort(vector<LoanRecipient> vec, int& reads, int& writes) {
    mergeSortTwoSortRec(vec, 0, vec.size() - 1, reads, writes);
}

void twoSort(vector<LoanRecipient> vec, int& reads, int& writes) {
    // This version of the first sorting algorithm sorts the vector by an identifier variable, then
    // calls the second sorting method, a version of bubble sort "bubbleSortTwoSort" to complete the
    // twoSort method
    selectionSortV1(vec, reads, writes);
}

void twoSortV2(vector<LoanRecipient> vec, int& reads, int& writes) {
    // This version of the first sorting algorithm sorts the vector by an identifier variable, then
    // calls the second sorting method, a version of merge sort "mergeSortTwoSort" to complete the
    // twoSortV2 method
    selectionSortV2(vec, reads, writes);
}

/** This sorting method was taken from Programming Algorithms at the following
 *  link: https://www.programmingalgorithms.com/algorithm/radix-sort/cpp/
 *  This method was altered to fit the LoanRecipient object type, as well as count
 *  the number of read and write actions occuring in the code. This link will also
 *  be present in the report. */

void RadixSort(vector<LoanRecipient> data, int count, int& reads, int& writes) {
    int i, j;
    vector<LoanRecipient> temp;
    for (int i = 0; i < data.size(); ++i) {
        temp.push_back(LoanRecipient()); // This vector loop will not be included in the number of reads or
                                         // writes for the counter, since this loop does not read from
                                         // the actual vector passed into the method.
    }

    for (int shift = 31; shift > -1; --shift)
    {
        j = 0;

        for (i = 0; i < count; ++i)
        {
            bool move = (data[i].getID() << shift) >= 0; // 1 read

            if (shift == 0 ? !move : move)
                data[i - j] = data[i]; // 1 read and 1 write
            else
                temp[j++] = data[i]; // 1 read and 1 write
            reads += 2;
            ++writes;
        }

        for (int i = 0; i < j; i++)
        {
            data[(count - j) + i] = temp[i]; // 1 read and 1 write
            ++reads;
            ++writes;
        }
    }
    temp.clear();
}

int main() {
        // Declared a vector of loan recipients
        vector<LoanRecipient> Loaners;
        string fileName = "../Training Data.csv/";
        // Initializing Datafiles from .csv file into LoanRecipient Vector Loaners
        if (getDataFromFile(fileName, Loaners)) {

            // This section of the code generates data from the sorting methods number of reads,
            // writes and times of completions as the number of items increases by a factor of 10.
            std::ofstream reads("readsSortingFile.csv");
            std::ofstream writes("writesSortingFile.csv");
            std::ofstream times("timeSortingFile.csv");
            reads << "Number Of Objects,Bubble Sort,Selection Sort,Insertion Sort,Merge Sort,Heap Sort,Two Sort,Two Sort V2,Radix Sort" << endl;
            writes << "Number Of Objects,Bubble Sort,Selection Sort,Insertion Sort,Merge Sort,Heap Sort,Two Sort,Two Sort V2,Radix Sort" << endl;
            times << "Number Of Objects,Bubble Sort,Selection Sort,Insertion Sort,Merge Sort,Heap Sort,Two Sort,Two Sort V2,Radix Sort" << endl;
            for (int r = 1; r < 11; ++r) {
                vector<LoanRecipient> loaners;
                for (int i = 0; i < r*100; ++i) {
                    loaners.push_back(LoanRecipient(Loaners[i].getID(), Loaners[i].getIncome(), Loaners[i].getAge(),
                                                    Loaners[i].getYearsOfExp(), Loaners[i].getMarStatus(),
                                                    Loaners[i].getHouseOwn(),Loaners[i].getCarOwn(),
                                                    Loaners[i].getJob(), Loaners[i].getCity(),
                                                    Loaners[i].getState(), Loaners[i].getCurrJobYrs(),
                                                    Loaners[i].getCurrHouseYrs(), Loaners[i].getRiskFlag()));
                }

                // In this vector here, I created a vector using a random number generator to
                // iterate values that may repeat, rather than just the iterator integer value of
                // the for-loop iteration. This random value is ensured to be within the given
                // range by dividing the value by 3 each time that it is greater than 1000. This
                // is done so that the number does not change favoring zero, along with decreasing
                // the integer by the least value possible. This is repeated until the vector is filled with
                // possibly-repeating integer values.

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsBubble = 0;
                int writesBubble = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startBubble;
                double bubbleTime;

                // The timer is started before either of the loops of the function begin
                startBubble = std::clock();

                bubbleSort(loaners, readsBubble, writesBubble);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                bubbleTime = ( std::clock() - startBubble ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsSelection = 0;
                int writesSelection = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startSelection;
                double selectionTime;

                // The timer is started before either of the loops of the function begin
                startSelection = std::clock();

                selectionSort(loaners, readsSelection, writesSelection);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                selectionTime = ( std::clock() - startSelection ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsInsertion = 0;
                int writesInsertion = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startInsertion;
                double insertionTime;

                // The timer is started before either of the loops of the function begin
                startInsertion = std::clock();

                insertionSort(loaners, readsInsertion, writesInsertion);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                insertionTime = ( std::clock() - startInsertion ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsMerge = 0;
                int writesMerge = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startMerge;
                double mergeTime;

                // The timer is started before either of the loops of the function begin
                startMerge = std::clock();

                mergeSort(loaners, readsMerge, writesMerge);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                mergeTime = ( std::clock() - startMerge ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsHeap = 0;
                int writesHeap = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startHeap;
                double heapTime;

                // The timer is started before either of the loops of the function begin
                startHeap = std::clock();

                heapSort(loaners, readsHeap, writesHeap);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                heapTime = ( std::clock() - startHeap ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsTwoSort = 0;
                int writesTwoSort = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startTwoSort;
                double twoSortTime;

                // The timer is started before either of the loops of the function begin
                startTwoSort = std::clock();

                twoSort(loaners, readsTwoSort, writesTwoSort);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                twoSortTime = ( std::clock() - startTwoSort ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsTwoSortV2 = 0;
                int writesTwoSortV2 = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startTwoSortV2;
                double twoSortV2Time;

                // The timer is started before either of the loops of the function begin
                startTwoSortV2 = std::clock();

                twoSortV2(loaners, readsTwoSortV2, writesTwoSortV2);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                twoSortV2Time = ( std::clock() - startTwoSortV2 ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsRadix = 0;
                int writesRadix = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startRadix;
                double radixTime;

                // The timer is started before either of the loops of the function begin
                startRadix = std::clock();

                RadixSort(loaners, r*100, readsRadix, writesRadix);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                radixTime = ( std::clock() - startRadix ) / (double) CLOCKS_PER_SEC;

                reads << r * 100 << "," << readsBubble << "," << readsSelection << "," << readsInsertion << "," << readsMerge <<
                "," << readsHeap << "," << readsTwoSort << "," << readsTwoSortV2 << "," << readsRadix << endl;
                writes << r * 100 << "," << writesBubble << "," << writesSelection << "," << writesInsertion << "," << writesMerge << "," <<
                writesHeap << "," << writesTwoSort << "," << writesTwoSortV2 << "," << writesRadix << endl;
                times << r * 100 << "," << bubbleTime << "," << selectionTime << "," << insertionTime << "," << mergeTime <<
                      "," << heapTime << "," << twoSortTime << "," << twoSortV2Time << "," << radixTime << endl;
            }

            // This section of the code generates data from the sorting methods number of reads,
            // writes and times of completions as the number of items increases by a factor of 2.
            std::ofstream readsTwoPow("readsSortingFileTwoPow.csv");
            std::ofstream writesTwoPow("writesSortingFileTwoPow.csv");
            std::ofstream timesTwoPow("timeSortingFileTwoPow.csv");
            readsTwoPow << "Number Of Objects,Bubble Sort,Selection Sort,Insertion Sort,Merge Sort,Heap Sort,Two Sort,Two Sort V2,Radix Sort" << endl;
            writesTwoPow << "Number Of Objects,Bubble Sort,Selection Sort,Insertion Sort,Merge Sort,Heap Sort,Two Sort,Two Sort V2,Radix Sort" << endl;
            timesTwoPow << "Number Of Objects,Bubble Sort,Selection Sort,Insertion Sort,Merge Sort,Heap Sort,Two Sort,Two Sort V2,Radix Sort" << endl;

            for (int r = 0; r < 6; ++r) {
                vector<LoanRecipient> loaners;
                int twoPOW = pow(2, r);
                for (int i = 0; i < twoPOW*100; ++i) {
                    loaners.push_back(LoanRecipient(Loaners[i].getID(), Loaners[i].getIncome(), Loaners[i].getAge(),
                                                    Loaners[i].getYearsOfExp(), Loaners[i].getMarStatus(),
                                                    Loaners[i].getHouseOwn(),Loaners[i].getCarOwn(),
                                                    Loaners[i].getJob(), Loaners[i].getCity(),
                                                    Loaners[i].getState(), Loaners[i].getCurrJobYrs(),
                                                    Loaners[i].getCurrHouseYrs(), Loaners[i].getRiskFlag()));
                }

                // In this vector here, I created a vector using a random number generator to
                // iterate values that may repeat, rather than just the iterator integer value of
                // the for-loop iteration. This random value is ensured to be within the given
                // range by dividing the value by 3 each time that it is greater than 1000. This
                // is done so that the number does not change favoring zero, along with decreasing
                // the integer by the least value possible. This is repeated until the vector is filled with
                // possibly-repeating integer values.

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsBubble = 0;
                int writesBubble = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startBubble;
                double bubbleTime;

                // The timer is started before either of the loops of the function begin
                startBubble = std::clock();

                bubbleSort(loaners, readsBubble, writesBubble);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                bubbleTime = ( std::clock() - startBubble ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsSelection = 0;
                int writesSelection = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startSelection;
                double selectionTime;

                // The timer is started before either of the loops of the function begin
                startSelection = std::clock();

                selectionSort(loaners, readsSelection, writesSelection);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                selectionTime = ( std::clock() - startSelection ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsInsertion = 0;
                int writesInsertion = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startInsertion;
                double insertionTime;

                // The timer is started before either of the loops of the function begin
                startInsertion = std::clock();

                insertionSort(loaners, readsInsertion, writesInsertion);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                insertionTime = ( std::clock() - startInsertion ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsMerge = 0;
                int writesMerge = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startMerge;
                double mergeTime;

                // The timer is started before either of the loops of the function begin
                startMerge = std::clock();

                mergeSort(loaners, readsMerge, writesMerge);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                mergeTime = ( std::clock() - startMerge ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsHeap = 0;
                int writesHeap = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startHeap;
                double heapTime;

                // The timer is started before either of the loops of the function begin
                startHeap = std::clock();

                heapSort(loaners, readsHeap, writesHeap);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                heapTime = ( std::clock() - startHeap ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsTwoSort = 0;
                int writesTwoSort = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startTwoSort;
                double twoSortTime;

                // The timer is started before either of the loops of the function begin
                startTwoSort = std::clock();

                twoSort(loaners, readsTwoSort, writesTwoSort);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                twoSortTime = ( std::clock() - startTwoSort ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsTwoSortV2 = 0;
                int writesTwoSortV2 = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startTwoSortV2;
                double twoSortV2Time;

                // The timer is started before either of the loops of the function begin
                startTwoSortV2 = std::clock();

                twoSortV2(loaners, readsTwoSortV2, writesTwoSortV2);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                twoSortV2Time = ( std::clock() - startTwoSortV2 ) / (double) CLOCKS_PER_SEC;

                std::shuffle(loaners.begin(), loaners.end(), std::default_random_engine());
                int readsRadix = 0;
                int writesRadix = 0;

                // For each timed Trial, I use a clock function clock_t to record the time it took for the tree to be created,
                // filled in, and each value searched in the tree. After this was done, the timer was stopped and the time of
                // both operations was printed to the user.
                clock_t startRadix;
                double radixTime;

                // The timer is started before either of the loops of the function begin
                startRadix = std::clock();

                RadixSort(loaners, twoPOW*100, readsRadix, writesRadix);

                // After each of the objects have been inserted and found in the tree, the timer is stopped and the
                // duration of time between start and stop is recorded.
                radixTime = ( std::clock() - startRadix ) / (double) CLOCKS_PER_SEC;

                readsTwoPow << twoPOW * 100 << "," << readsBubble << "," << readsSelection << "," << readsInsertion << "," << readsMerge <<
                      "," << readsHeap << "," << readsTwoSort << "," << readsTwoSortV2 << "," << readsRadix << endl;
                writesTwoPow << twoPOW * 100 << "," << writesBubble << "," << writesSelection << "," << writesInsertion << "," << writesMerge << "," <<
                       writesHeap << "," << writesTwoSort << "," << writesTwoSortV2 << "," << writesRadix << endl;
                timesTwoPow << twoPOW * 100 << "," << bubbleTime << "," << selectionTime << "," << insertionTime << "," << mergeTime <<
                      "," << heapTime << "," << twoSortTime << "," << twoSortV2Time << "," << radixTime << endl;
            }

// It should be noted that each of these shuffles are randomized. Because of this, each subsequent running of the
// program could cause different number of reads and writes from each program, creating a possible source of
// error while also mitigating the chance of one sorting method benefiting from the initial shuffling of the
// vector. While the number of reads and writes should increase by a consistent factor, if the number of reads
// and writes are small enough this factor could possibly change the results of our tests dramatically. After
// initial testing results were returned however, it appears that the number of results returned mitigate the
// possibility of the general trend being altered by this.

            return 0;
        }
    }
