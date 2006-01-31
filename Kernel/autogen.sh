#!/bin/sh -x

# Print every line before execution
set -e

# First cache the command names in variables. If you want to
# override the names, simply set the variables before calling this
# script.
: ${LIBTOOLIZE=libtoolize}
: ${ACLOCAL=aclocal}
: ${AUTOHEADER=autoheader}
: ${AUTOMAKE=automake}
: ${AUTOCONF=autoconf}

# Call the commands one by one
${LIBTOOLIZE} -f --automake
${ACLOCAL} -I macros ${ACLOCAL_FLAGS}
${AUTOHEADER}
${AUTOMAKE} --add-missing
${AUTOCONF}

# Successfully finished.
echo "Now you can run ./configure"
