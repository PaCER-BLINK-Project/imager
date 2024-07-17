#ifndef __BLINK_PIPELINE_UTILS__
#define __BLINK_PIPELINE_UTILS__

#include <string>

namespace blink::imager {
    /**
     * @brief Checks whether a path corresponds to an existing directory.
     */
    bool dir_exists(const std::string& path);



    /**
     * @brief Creates a directory at the given path.
     *
     * If a directory already exists, nothing is done.
     */
    void create_directory(const std::string& path);

}
#endif