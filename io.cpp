#include "io.h"

#include <iostream>
#include <set>
#include <algorithm>
#include <cstdlib>

// Remove all characters 'c' from string 'str'.
// Also output the number of characters removed
int remove_all_chars_and_count(
  char* str, 
  const char c) 

{
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

// Check if a given string is a number
bool legal_int(
  const char *str) 

{
  while (*str)
    if (!isdigit(*str++))
      return false;
  return true;
}