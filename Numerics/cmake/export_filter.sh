gawk '$2 ~ /(D|T|B)/ { print $3 }' | sed 1iEXPORTS
