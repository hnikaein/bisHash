#include "time_profile.h"
#include <chrono>
#include <vector>
#include <map>
#include <pthread.h>


using namespace std;

auto times = map<string, vector<chrono::milliseconds>>();
pthread_mutex_t times_mutex;

void add_time_c(const string &caller_name) {
    pthread_mutex_lock(&times_mutex);
    if (!times.count(caller_name))
        times[caller_name] = vector<chrono::milliseconds>();
    times[caller_name].push_back(duration_cast<chrono::milliseconds>(chrono::system_clock::now().time_since_epoch()));
    pthread_mutex_unlock(&times_mutex);
}

string get_times_str_c(const string &caller_name, bool free_space) {
    string res;
    pthread_mutex_lock(&times_mutex);
    for (int i = 1; i < times[caller_name].size(); i++)
        res += to_string((times[caller_name][i] - times[caller_name][i - 1]).count()) + " ";
    pthread_mutex_unlock(&times_mutex);
    if (free_space)
        erase_times_c(caller_name);
    return res;
}

long last_time_c(const string &caller_name) {
    pthread_mutex_lock(&times_mutex);
    unsigned long siz = times[caller_name].size();
    long result = times[caller_name][siz - 1].count() - times[caller_name][siz - 2].count();
    pthread_mutex_unlock(&times_mutex);
    return result;
}

void erase_times_c(const string &caller_name) {
    pthread_mutex_lock(&times_mutex);
    times.erase(caller_name);
    pthread_mutex_unlock(&times_mutex);
}