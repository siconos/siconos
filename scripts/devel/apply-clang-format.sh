#!/bin/bash

# Check all of the arguments first to make sure they're all directories
for dir in "$@"; do
    if [ ! -d "${dir}" ]; then
        echo "${dir} is not a directory"
        usage
    fi
done

# Find a dominating file, starting from a given directory and going up.
find-dominating-file() {
    if [ -r "$1"/"$2" ]; then
        return 0
    fi
    if [ "$1" = "/" ]; then
        return 1
    fi
    find-dominating-file "$(realpath "$1"/..)" "$2"
    return $?
}

# Run clang-format -i on all of the things
for dir in "$@"; do
    pushd "${dir}" &>/dev/null
    if ! find-dominating-file . .clang-format; then
        echo "Failed to find dominating .clang-format starting at $PWD"
        continue
    fi
    find . \
         \( -name '*.c' \
         -o -name '*.cpp' \
         -o -name '*.h' \
         -o -name '*.hpp' \) \
         -exec clang-format -i '{}' \;
    #popd &>/dev/null
done