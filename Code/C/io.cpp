#include "io.h"

int remove_all_chars_and_count(char* str, char c) {
    int removed = 0;
    char *pr = str, *pw = str;
    while (*pr) {
        *pw = *pr++;
        removed += (*pw == c);
        pw += (*pw != c);
    }
    *pw = '\0';
    return removed;
}

// function to check is a char * is an integer
bool legal_int(char *str) {
    while (*str)
        if (!isdigit(*str++))
            return false;
    return true;
}