#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <string>
#include <string.h>

namespace blink {

    namespace imager {

        void create_directory(const std::string& path){
            char tmp[1024];
            char *p = NULL;
            size_t len;

            snprintf(tmp, sizeof(tmp),"%s", path.c_str());
            len = strlen(tmp);
            if (tmp[len - 1] == '/')
                tmp[len - 1] = 0;
            for (p = tmp + 1; *p; p++)
                if (*p == '/') {
                    *p = 0;
                    mkdir(tmp, S_IRWXU);
                    *p = '/';
                }
            mkdir(tmp, S_IRWXU);
        }

        bool dir_exists(const std::string& path){
            DIR* dir = opendir(path.c_str());
            if (dir) {
                closedir(dir);
                return true;
            } else {
                return false;
            }
        }
    }
}