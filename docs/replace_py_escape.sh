#!/bin/bash


sed -e 's/\\\\/\\n/g; s/\\r/\\\\r/g; s/\\a/\\\\a/g;s/\\t/\\\\t/g; s/\\v/\\\\v/g; s/\\x/\\\\x/g; s/\\f/\\\\f/g; s/\\b/\\\\b/g' $1  > $2
# sed 's/\\x/\\\\x/g' $2  > $2
# sed 's/\\v/\\\\v/g' $2  > $2
