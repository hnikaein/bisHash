#include <vector>

#ifndef MULTIPROC_H
#define MULTIPROC_H


int multiproc(int thread_count, int (*start_routine)(const int), int end, int from = 0);

int multiproc(int thread_count, int (*start_routine)(const int), const std::vector<int> &ids);

#endif //MULTIPROC_H
