#include <string>
#include <mutex>

#ifndef LOGGER_H
#define LOGGER_H

#define FORMAT_LENGTH 1000

class Logger {
public:
    enum [[maybe_unused]] LogLevel {
        ALL = 1000, DEBUGL2 = 6, DEBUG = 5, INFO = 4, WARN = 3, ERROR = 2, FATAL = 1, OFF = 0
    };

    explicit Logger(int log_level);

    void debugl2_noheader(const char *format, ...);

    void debugl2(const char *format, ...);

    void debug(const char *format, ...);

    void info(const char *format, ...);

    [[maybe_unused]] void warn(const char *format, ...);

    [[maybe_unused]] void error(const char *format, ...);

    [[maybe_unused]] void fatal(const char *format, ...);

    void log(const std::string &s, LogLevel log_level, bool noheader = false);

    static std::string formatString(const char *format, ...);

    LogLevel log_level = OFF;
private:
    std::mutex mtx;

    static std::string formatString(const char *format, va_list args);
};

#endif //LOGGER_H