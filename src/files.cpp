#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>

namespace blink {

    bool dir_exists(const std::string& path){
        DIR* dir = opendir(path.c_str());
        if (dir) {
            closedir(dir);
            return true;
        } else {
            return false;
        }
    }


    bool create_directory(const std::string& path){
        if(!dir_exists(path.c_str())){
            if(mkdir(path.c_str(), 0775) != 0) return false;
            return true;
        }
        return false;
    }
}