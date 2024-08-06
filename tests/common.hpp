#ifndef __COMMON_H__
#define __COMMON_H__

#include <exception>
#include <string>

#define ENV_DATA_ROOT_DIR "BLINK_TEST_DATADIR"

class TestFailed : public std::exception {
    private:
    std::string message;

    public:
    TestFailed(const char *msg)  noexcept : std::exception(),  message {msg} {}
    virtual ~TestFailed() {};
    virtual const char* what() const noexcept override {return message.c_str();}
};


#endif
